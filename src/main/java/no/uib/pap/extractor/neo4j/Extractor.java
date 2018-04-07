package no.uib.pap.extractor.neo4j;

import com.google.common.collect.*;
import no.uib.pap.model.*;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.neo4j.driver.v1.Record;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static no.uib.pap.model.Error.ERROR_READING_VEP_TABLES;
import static no.uib.pap.model.Error.sendError;

public class Extractor {

    // Separate sets of objects are created to have more than one attribute of the
    // objects like, imapReactions or iPathways
    // Other objects like genes, proteins don't need a separate list, because the
    // identifier is the only attribute used.

    private static ImmutableMap<String, Reaction> imapReactions; // Reaction stId to Reaction displayName
    private static ImmutableMap<String, Pathway> iPathways; // Pathway stId to Pathway instance
    private static ImmutableMap<String, String> iProteins; // Protein accession (UniProt) to name
    private static ImmutableSetMultimap<String, String> imapGenesToProteins = null;
    private static ImmutableSetMultimap<String, String> imapEnsemblToProteins = null;
    private static ImmutableSetMultimap<String, String> imapPhysicalEntitiesToReactions = null;
    private static ImmutableSetMultimap<String, Proteoform> imapPhysicalEntitiesToProteoforms = null;
    private static ImmutableSetMultimap<Proteoform, String> imapProteoformsToPhysicalEntities = null;
    private static ImmutableSetMultimap<String, String> imapProteinsToReactions = null;
    private static ImmutableSetMultimap<String, String> imapReactionsToPathways = null;
    private static ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways = null;
    private static ImmutableSetMultimap<String, Proteoform> imapProteinsToProteoforms = null;
    private static ImmutableSetMultimap<Proteoform, String> imapProteoformsToReactions = null;
    private static ImmutableSetMultimap<String, String> imapProteinsToComplexes = null;
    private static ImmutableSetMultimap<String, String> imapComplexesToProteins = null;
    private static ImmutableSetMultimap<String, String> imapSetsToProteins = null;
    private static ImmutableSetMultimap<String, String> imapProteinsToSets = null;
    private static ImmutableSetMultimap<Proteoform, String> imapProteoformsToComplexes = null;
    private static ImmutableSetMultimap<String, Proteoform> imapComplexesToProteoforms = null;
    private static ImmutableSetMultimap<Proteoform, String> imapProteoformsToSets = null;
    private static ImmutableSetMultimap<String, Proteoform> imapSetsToProteoforms = null;

    private static ImmutableSetMultimap<String, String> imapRsIdsToProteins = null; // An array of multimaps, one for each chromosome. Added one extra to use natural 1-based numbering.
    private static ImmutableSetMultimap<Long, String> imapChrBpToProteins = null;

    private static final int rsidColumnIndex = 2;
    // private static final int ensemblColumnIndex = 4;
    private static final int swissprotColumnIndex = 5;
    // private static final int nearestGeneColumnIndex = 7;

    static String outputPath = "resources/";

    public static void main(String[] args) {

        System.out.println("The working directory is: " + System.getProperty("user.dir"));

        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        imapProteinsToReactions = getProteinsToReactions();
        System.out.println("Finished map proteins to imapReactions.");

        for (int chr = 1; chr <= 22; chr++) {
            imapRsIdsToProteins = getRsIdsToProteins(chr);
            System.out.println("Finished map rsids to proteins, chromosome " + chr);
        }

        for (int chr = 1; chr <= 22; chr++) {
            imapChrBpToProteins = getChrBpToProteins(chr);
            System.out.println("Finished map chrBp to proteins, chromosome " + chr);
        }

        imapGenesToProteins = getGenesToProteins();
        System.out.println("Finished map genes to proteins.");

        imapEnsemblToProteins = getEnsemblToProteins();
        System.out.println("Finished map ensembl to proteins.");

        imapProteoformsToReactions = getProteoformsToReactions();
        System.out.println("Finished map proteoforms to reactions.");

        imapReactionsToPathways = getReactonsToPathways();
        System.out.println("Finished map imapReactions to iPathways.");

        imapPathwaysToTopLevelPathways = getPathwaysToTopLevelPathways();
        System.out.println("Finished map iPathways to top level iPathways.");

        iProteins = getProteinNames();
        System.out.println("Finished getting the protein names.");

        imapProteinsToComplexes = getProteinsToComplexes();
        System.out.println("Finished map proteins to complexes.");

        imapComplexesToProteins = getComplexesToProteins();
        System.out.println("Finished map of complexes to proteins.");

        imapProteoformsToComplexes = getProteoformsToComplexes();
        System.out.println("Finished map of proteoforms to complexes");

        imapComplexesToProteoforms = getComplexesToProteoforms();
        System.out.println("Finished map of complexes to proteoforms");

        imapSetsToProteins = getSetsToProteins();
        System.out.println("Finished map of sets to proteins.");

        imapProteinsToSets = getProteinsToSets();
        System.out.println("Finished map proteins to sets");

        imapSetsToProteoforms = getSetsToProteoforms();
        System.out.println("Finished map sets to proteoforms");

        imapProteoformsToSets = getProteoformsToSets();
        System.out.println("Finished map proteoforms to sets");
    }

    private static void getComplexComponents() {

        if (imapPhysicalEntitiesToProteoforms == null) {
            getPhysicalEntitiesToProteoforms();
        }

        ImmutableSetMultimap.Builder<String, String> builderComplexesToProteins = new ImmutableSetMultimap.Builder<>();
        ImmutableSetMultimap.Builder<String, String> builderProteinsToComplexes = new ImmutableSetMultimap.Builder<>();
        ImmutableSetMultimap.Builder<String, Proteoform> builderComplexesToProteoforms = new ImmutableSetMultimap.Builder<>();
        ImmutableSetMultimap.Builder<Proteoform, String> builderProteoformsToComplexes = new ImmutableSetMultimap.Builder<>();

        //Query the database to fill the data structure
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_COMPLEX_COMPONENTS);
        for (Record record : resultList) {
            builderComplexesToProteins.put(record.get("complex").asString(), record.get("protein").asString());
            builderProteinsToComplexes.put(record.get("protein").asString(), record.get("complex").asString());

            for (Proteoform proteoform : imapPhysicalEntitiesToProteoforms.get(record.get("physicalEntity").asString())) {
                builderComplexesToProteoforms.put(record.get("complex").asString(), proteoform);
                builderProteoformsToComplexes.put(proteoform, record.get("complex").asString());
            }
        }

        imapComplexesToProteins = builderComplexesToProteins.build();
        imapProteinsToComplexes = builderProteinsToComplexes.build();
        imapComplexesToProteoforms = builderComplexesToProteoforms.build();
        imapProteoformsToComplexes = builderProteoformsToComplexes.build();

        storeSerialized(imapComplexesToProteins, outputPath + "imapComplexesToProteins.gz");
        storeSerialized(imapProteinsToComplexes, outputPath + "imapProteinsToComplexes.gz");
        storeSerialized(imapComplexesToProteoforms, outputPath + "imapComplexesToProteoforms.gz");
        storeSerialized(imapProteoformsToComplexes, outputPath + "imapProteoformsToComplexes.gz");
    }

    public static ImmutableSetMultimap<String, String> getComplexesToProteins() {
        if (imapComplexesToProteins == null) {
            getComplexComponents();
        }

        return imapComplexesToProteins;
    }


    public static ImmutableSetMultimap<String, String> getProteinsToComplexes() {
        if (imapProteinsToComplexes == null) {
            getComplexComponents();
        }
        return imapProteinsToComplexes;
    }

    public static ImmutableSetMultimap<String, Proteoform> getComplexesToProteoforms() {
        if (imapComplexesToProteoforms == null) {
            getComplexComponents();
        }

        return imapComplexesToProteoforms;
    }

    public static ImmutableSetMultimap<Proteoform, String> getProteoformsToComplexes() {
        if (imapProteoformsToComplexes == null) {
            getComplexComponents();
        }
        return imapProteoformsToComplexes;
    }

    public static ImmutableSetMultimap<String, String> getSetMembersAndCandidates() {
        if(imapPhysicalEntitiesToProteoforms == null){
            getPhysicalEntitiesToProteoforms();
        }

        ImmutableSetMultimap.Builder<String, String> builderSetsToProteins = new ImmutableSetMultimap.Builder<>();
        ImmutableSetMultimap.Builder<String, String> builderProteinsToSets = new ImmutableSetMultimap.Builder<>();
        ImmutableSetMultimap.Builder<String, Proteoform> builderSetsToProteoforms = new ImmutableSetMultimap.Builder<>();
        ImmutableSetMultimap.Builder<Proteoform, String> builderProteoformsToSets = new ImmutableSetMultimap.Builder<>();

        //Query the database to fill the data structure
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_SET_MEMBERS_AND_CANDIDATES);
        for (Record record : resultList) {
            builderSetsToProteins.put(record.get("set").asString(), record.get("protein").asString());
            builderProteinsToSets.put(record.get("protein").asString(), record.get("set").asString());
            for (Proteoform proteoform : imapPhysicalEntitiesToProteoforms.get(record.get("physicalEntity").asString())) {
                builderSetsToProteoforms.put(record.get("set").asString(), proteoform);
                builderProteoformsToSets.put(proteoform, record.get("set").asString());
            }
        }

        imapSetsToProteins = builderSetsToProteins.build();
        imapProteinsToSets = builderProteinsToSets.build();
        imapSetsToProteoforms = builderSetsToProteoforms.build();
        imapProteoformsToSets = builderProteoformsToSets.build();

        storeSerialized(imapSetsToProteins, outputPath + "imapSetsToProteins.gz");
        storeSerialized(imapProteinsToSets, outputPath + "imapProteinsToSets.gz");
        storeSerialized(imapSetsToProteoforms, outputPath + "imapSetsToProteoforms.gz");
        storeSerialized(imapProteoformsToSets, outputPath + "imapProteoformsToSets.gz");

        return imapSetsToProteins;
    }

    public static ImmutableSetMultimap<String, String> getSetsToProteins() {
        if (imapSetsToProteins == null) {
            getSetMembersAndCandidates();
        }

        return imapSetsToProteins;
    }


    public static ImmutableSetMultimap<String, String> getProteinsToSets() {
        if (imapProteinsToSets == null) {
            getSetMembersAndCandidates();
        }
        return imapProteinsToSets;
    }

    public static ImmutableSetMultimap<String, Proteoform> getSetsToProteoforms() {
        if (imapSetsToProteoforms == null) {
            getSetMembersAndCandidates();
        }

        return imapSetsToProteoforms;
    }

    public static ImmutableSetMultimap<Proteoform, String> getProteoformsToSets() {
        if (imapProteoformsToSets == null) {
            getSetMembersAndCandidates();
        }
        return imapProteoformsToSets;
    }


    public static ImmutableSetMultimap<String, String> getRsIdsToProteins(int chr) {

        if (imapProteinsToReactions == null) {
            getProteinsToReactions();
        }

        ImmutableSetMultimap.Builder<String, String> builderRsIdsToProteins = ImmutableSetMultimap.<String, String>builder();

        // Traverse all the vepTables

        System.out.println("Scanning vepTable for chromosome " + chr);
        try {
            BufferedReader br = getBufferedReaderFromResource(chr + ".gz");
            br.readLine(); // Read header line

            for (String line; (line = br.readLine()) != null; ) {

                Multimap<Snp, String> snpToSwissprotMap = getSNPAndSwissProtFromVep(line);

                for (Map.Entry<Snp, String> snpToSwissprotPair : snpToSwissprotMap.entries()) {
                    Snp snp = snpToSwissprotPair.getKey();
                    String protein = snpToSwissprotPair.getValue();

                    if (!protein.equals("NA")) {
                        if (imapProteinsToReactions.containsKey(protein)) {
                            builderRsIdsToProteins.put(snp.getRsid(), protein);
                        }
                    }
                }
            }
        } catch (IOException ex) {
            sendError(ERROR_READING_VEP_TABLES, chr);
        }


        imapRsIdsToProteins = builderRsIdsToProteins.build();

        storeSerialized(imapRsIdsToProteins, outputPath + "imapRsIdsToProteins" + chr + ".gz");

        return imapRsIdsToProteins;
    }

    public static ImmutableSetMultimap<Long, String> getChrBpToProteins(int chr) {

        if (imapProteinsToReactions == null) {
            getProteinsToReactions();
        }

        ImmutableSetMultimap.Builder<Long, String> builderChrBpToProteins = ImmutableSetMultimap.<Long, String>builder();

        // Traverse all the vepTables

        System.out.println("Scanning vepTable for chromosome " + chr);
        try {
            BufferedReader br = getBufferedReaderFromResource(chr + ".gz");
            br.readLine(); // Read header line

            for (String line; (line = br.readLine()) != null; ) {

                Multimap<Snp, String> snpToSwissprotMap = getSNPAndSwissProtFromVep(line);

                for (Map.Entry<Snp, String> snpToSwissprotPair : snpToSwissprotMap.entries()) {
                    Snp snp = snpToSwissprotPair.getKey();
                    String protein = snpToSwissprotPair.getValue();

                    if (!protein.equals("NA")) {
                        if (imapProteinsToReactions.containsKey(protein)) {
                            builderChrBpToProteins.put(snp.getBp(), protein);
                        }
                    }
                }
            }
        } catch (IOException ex) {
            sendError(ERROR_READING_VEP_TABLES, chr);
        }

        imapChrBpToProteins = builderChrBpToProteins.build();

        storeSerialized(imapChrBpToProteins, outputPath + "imapChrBpToProteins" + chr + ".gz");
        return imapChrBpToProteins;
    }

    public static ImmutableSetMultimap<String, Proteoform> getProteinsToProteoforms() {

        // Make sure the full list of proteoforms is loaded
        if (imapProteoformsToPhysicalEntities == null) {
            getPhysicalEntitiesToProteoforms();
        }

        ImmutableSetMultimap.Builder<String, Proteoform> builderProteinsToProteoforms = ImmutableSetMultimap
                .<String, Proteoform>builder();

        //Traverse the list of proteoforms to get the protein accessions
        for (Proteoform proteoform : imapProteoformsToPhysicalEntities.keySet()) {
            builderProteinsToProteoforms.put(proteoform.getUniProtAcc(), proteoform);
        }

        imapProteinsToProteoforms = builderProteinsToProteoforms.build();
        storeSerialized(imapProteinsToProteoforms, outputPath + "imapProteinsToProteoforms.gz");

        return imapProteinsToProteoforms;
    }

    public static ImmutableSetMultimap<String, Proteoform> getPhysicalEntitiesToProteoforms() {

        ImmutableSetMultimap.Builder<String, Proteoform> builderPhysicalEntitiesToProteoforms = ImmutableSetMultimap.<String, Proteoform>builder();
        ImmutableSetMultimap.Builder<Proteoform, String> builderProteoformsToPhysicalEntities = ImmutableSetMultimap.<Proteoform, String>builder();

        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PROTEOFORMS_TO_PHYSICALENTITIES);
        for (Record record : resultList) {

            // Create the proteoform instance using the isoform and ptms.
            String isoform = record.get("isoform").asString();
            List<Pair<String, Long>> ptms = new ArrayList<>();

            for (Object ptm : record.get("ptms").asList()) {
                String[] parts = ptm.toString().split(":");
                String mod = parts[0];
                Long site = Proteoform.interpretCoordinateFromStringToLong(parts[1]);
                ptms.add(new MutablePair<>(mod, site));
            }

            Proteoform proteoform = new Proteoform(isoform, ptms);

            for (Object physicalEntity : record.get("peSet").asList()) {
                builderPhysicalEntitiesToProteoforms.put(physicalEntity.toString(), proteoform);
                builderProteoformsToPhysicalEntities.putAll(proteoform, physicalEntity.toString());
            }
        }

        imapPhysicalEntitiesToProteoforms = builderPhysicalEntitiesToProteoforms.build();
        imapProteoformsToPhysicalEntities = builderProteoformsToPhysicalEntities.build();

        return imapPhysicalEntitiesToProteoforms;
    }

    public static ImmutableSetMultimap<Proteoform, String> getProteoformsToPhysicalEntities() {
        if (imapProteoformsToPhysicalEntities == null) {
            getPhysicalEntitiesToProteoforms();
        }
        return imapProteoformsToPhysicalEntities;
    }

    /**
     * Get list of imapReactions
     */
    public static ImmutableMap<String, Reaction> getReactions() {

        if (imapPhysicalEntitiesToProteoforms == null) {
            getPhysicalEntitiesToProteoforms();
        }

        ImmutableMap.Builder<String, Reaction> builderReactions = ImmutableMap.<String, Reaction>builder();

        // Fill the reactions stId and displayName
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_ALL_REACTIONS);
        for (Record record : resultList) {
            Reaction reaction = new Reaction(record.get("stId").asString(), record.get("displayName").asString());
            builderReactions.put(record.get("stId").asString(), reaction);
        }

        imapReactions = builderReactions.build();

        // Fill the reaction participants
        resultList = ConnectionNeo4j.query(ReactomeQueries.GET_REACTION_PARTICIPANTS);
        for (Record record : resultList) {
            String role = record.get("role").asString().toUpperCase();
            if (!imapReactions.containsKey(record.get("reaction").asString())) {
                System.out.println("Missing reaction: " + record.get("reaction").asString());
            }
            imapReactions.get(record.get("reaction").asString()).addParticipant(record.get("protein").asString(), Role.valueOf(role));
            for (Proteoform proteoform : imapPhysicalEntitiesToProteoforms.get(record.get("physicalEntity").asString())) {
                imapReactions.get(record.get("reaction").asString()).addParticipant(proteoform, Role.valueOf(role));
            }
        }

        // Serialize list of imapReactions
        storeSerialized(imapReactions, outputPath + "imapReactions.gz");
        return imapReactions;
    }

    /**
     * Get list of iPathways
     */
    public static ImmutableMap<String, Pathway> getPathways() {

        ImmutableMap.Builder<String, Pathway> builderPathways = ImmutableMap.<String, Pathway>builder();

        // Query the database and fill the data structure
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_ALL_PATHWAYS);
        for (Record record : resultList) {

            // Get first the attributes that can be obtained directly
            Pathway pathway = new Pathway(record.get("stId").asString(), record.get("displayName").asString());
            pathway.setNumEntitiesTotal(record.get("numEntitiesTotal").asInt());
            pathway.setNumReactionsTotal(record.get("numReactionsTotal").asInt());
            pathway.setNumReactionsTotal(record.get("numProteoformsTotal").asInt());
            builderPathways.put(pathway.getStId(), pathway);
        }

        iPathways = builderPathways.build();

        // Serialize iPathways
        storeSerialized(iPathways, outputPath + "iPathways.gz");
        return iPathways;
    }

    public static ImmutableSetMultimap<String, String> getGenesToProteins() {

        // Query the database and fill the data structure
        ImmutableSetMultimap.Builder<String, String> builderGenesToProteins = ImmutableSetMultimap
                .<String, String>builder();
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_GENES_TO_PROTEINS);
        for (Record record : resultList) {
            builderGenesToProteins.put(record.get("gene").asString(), record.get("protein").asString());
        }
        imapGenesToProteins = builderGenesToProteins.build();

        storeSerialized(imapGenesToProteins, outputPath + "imapGenesToProteins.gz");
        return imapGenesToProteins;
    }

    public static ImmutableSetMultimap<String, String> getEnsemblToProteins() {

        // Query the database and fill the data structure
        ImmutableSetMultimap.Builder<String, String> builderEnsemblToProteins = ImmutableSetMultimap
                .<String, String>builder();
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_ENSEMBL_TO_PROTEINS);
        for (Record record : resultList) {
            builderEnsemblToProteins.put(record.get("ensembl").asString(), record.get("uniprot").asString());
        }
        imapEnsemblToProteins = builderEnsemblToProteins.build();

        storeSerialized(imapEnsemblToProteins, outputPath + "imapEnsemblToProteins.gz");
        return imapEnsemblToProteins;
    }

    /**
     * Get physical entities to reactions mapping
     */
    public static ImmutableSetMultimap<String, String> getPhysicalEntitiesToReactions() {

        ImmutableSetMultimap.Builder<String, String> builderPhysicalEntitiesToReactions = ImmutableSetMultimap
                .<String, String>builder();

        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PHYSICALENTITIES_TO_REACTIONS);
        for (Record record : resultList) {
            builderPhysicalEntitiesToReactions.put(
                    record.get("physicalEntity").asString(),
                    record.get("reaction").asString()
            );
        }

        imapPhysicalEntitiesToReactions = builderPhysicalEntitiesToReactions.build();
        return imapPhysicalEntitiesToReactions;
    }


    /**
     * Creates mapping from proteins to proteoforms and mapping from proteoforms to reactions.
     * Note that the proteins are identified by uniprot accession without the isoform, but the proteoforms include the isoform.
     *
     * @return
     */
    public static ImmutableSetMultimap<Proteoform, String> getProteoformsToReactions() {
        if (imapPhysicalEntitiesToReactions == null) {
            getPhysicalEntitiesToReactions();
        }
        if (imapProteoformsToPhysicalEntities == null) {
            getPhysicalEntitiesToProteoforms();
        }

        ImmutableSetMultimap.Builder<Proteoform, String> builderProteoformsToReactions = ImmutableSetMultimap
                .<Proteoform, String>builder();

        for (Map.Entry<Proteoform, String> proteoformToPhysicalEntity : imapProteoformsToPhysicalEntities.entries()) {

            for (String reaction : imapPhysicalEntitiesToReactions.get(proteoformToPhysicalEntity.getValue())) {
                builderProteoformsToReactions.put(proteoformToPhysicalEntity.getKey(), reaction);
            }
        }

        imapProteoformsToReactions = builderProteoformsToReactions.build();

        storeSerialized(imapProteoformsToReactions, outputPath + "imapProteoformsToReactions.gz");

        return imapProteoformsToReactions;
    }

    /**
     * Get mapping from proteins to imapReactions
     */
    public static ImmutableSetMultimap<String, String> getProteinsToReactions() {

        if (imapReactions == null) {
            getReactions();
        }

        // Query the database and fill the data structure
        ImmutableSetMultimap.Builder<String, String> builderProteinsToReactions = ImmutableSetMultimap
                .<String, String>builder();
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PROTEINS_TO_REACTIONS);

        for (Record record : resultList) {
            if (record.get("reaction").asString().equals("null")) {
                builderProteinsToReactions.putAll(record.get("protein").asString(), new ArrayList<String>());
            } else {
                builderProteinsToReactions.put(record.get("protein").asString(), record.get("reaction").asString());
            }
        }

        imapProteinsToReactions = builderProteinsToReactions.build();
        storeSerialized(imapProteinsToReactions, outputPath + "imapProteinsToReactions.gz");
        return imapProteinsToReactions;
    }

    public static ImmutableMap<String, String> getProteinNames() {
        ImmutableMap.Builder<String, String> proteinsBuilder = ImmutableMap.<String, String>builder();

        String name = "";
        String query = ReactomeQueries.GET_PROTEIN_NAMES;
        List<Record> resultList = ConnectionNeo4j.query(query);

        for (Record record : resultList) {
            for (Object entry : record.get("description").asList()) {
//                System.out.println("----" + entry.toString() + "---");
                name = record.get("description").asList().toString();
                String[] parts = name.split("(recommendedName:)|(alternativeName:)|(component recommendedName:)");
                if (parts.length <= 1) {
                    name = record.get("displayName").asString();
                } else {
                    name = parts[1].trim();
                }
                proteinsBuilder.put(record.get("identifier").asString(), name);
//                System.out.println(name + "\n");
                break;
            }
        }

        iProteins = proteinsBuilder.build();
        storeSerialized(iProteins, outputPath + "iProteins.gz");
        return iProteins;
    }

    /**
     * Get mapping from imapReactions to iPathways
     */
    private static ImmutableSetMultimap<String, String> getReactonsToPathways() {

        if (imapReactions == null) {
            getReactions();
        }
        if (iPathways == null) {
            getPathways();
        }

        // Query the database and fill the data structure
        ImmutableSetMultimap.Builder<String, String> builderReactionsToPathways = ImmutableSetMultimap
                .<String, String>builder();
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_REACTIONS_TO_PATHWAYS);
        for (Record record : resultList) {
            builderReactionsToPathways.put(record.get("reaction").asString(), record.get("pathway").asString());
        }

        imapReactionsToPathways = builderReactionsToPathways.build();

        storeSerialized(imapReactionsToPathways, outputPath + "imapReactionsToPathways.gz");
        return imapReactionsToPathways;
    }

    /**
     * Get mapping from iPathways to top level iPathways
     */
    public static ImmutableSetMultimap<String, String> getPathwaysToTopLevelPathways() {

        if (iPathways.size() == 0) {
            getPathways();
        }

        // Query the database and fill the data structure
        ImmutableSetMultimap.Builder<String, String> builderPathwaysToTopLevelPathways = ImmutableSetMultimap
                .<String, String>builder();
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GEP_MAP_PATHWAYS_TO_TOPLEVELPATHWAYS);
        for (Record record : resultList) {
            builderPathwaysToTopLevelPathways.put(record.get("pathway").asString(),
                    record.get("topLevelPathway").asString());
        }

        imapPathwaysToTopLevelPathways = builderPathwaysToTopLevelPathways.build();

        storeSerialized(imapPathwaysToTopLevelPathways, outputPath + "imapPathwaysToTopLevelPathways.gz");
        return imapPathwaysToTopLevelPathways;
    }

    public static Multimap<Snp, String> getSNPAndSwissProtFromVep(String line) {
        ImmutableSetMultimap.Builder<Snp, String> builder = ImmutableSetMultimap.<Snp, String>builder();
        String[] fields = line.split(" ");
        Integer chr = Integer.valueOf(fields[0]);
        Long bp = Long.valueOf(fields[1]);

        String[] rsids = fields[rsidColumnIndex].split(",");
        String[] uniprots = fields[swissprotColumnIndex].split(",");

        for (String rsid : rsids) {
            for (String uniprot : uniprots) {
                Snp snp = new Snp(chr, bp, rsid);
                builder.put(snp, uniprot);
            }
        }
        return builder.build();
    }

    static BufferedReader getBufferedReaderFromResource(String fileName) throws FileNotFoundException, IOException {

        BufferedReader br = null;
        InputStream fileStream = ClassLoader.getSystemResourceAsStream(fileName);
        InputStream gzipStream = new GZIPInputStream(fileStream);
        Reader decoder = new InputStreamReader(gzipStream);
        br = new BufferedReader(decoder);

        return br;
    }

    private static void storeSerialized(Object obj, String fileName) {
        FileOutputStream fos = null;
        GZIPOutputStream gz;
        ObjectOutputStream oos;
        try {
            fos = new FileOutputStream(fileName);
            gz = new GZIPOutputStream(fos);
            oos = new ObjectOutputStream(gz);
            oos.writeObject(obj);
            oos.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
