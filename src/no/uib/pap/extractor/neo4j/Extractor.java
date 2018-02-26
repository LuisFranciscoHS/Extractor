package no.uib.pap.extractor.neo4j;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Multimap;
import no.uib.pap.model.Pathway;
import no.uib.pap.model.Proteoform;
import no.uib.pap.model.Snp;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.neo4j.driver.internal.value.StringValue;
import org.neo4j.driver.v1.Record;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static no.uib.pap.model.Error.ERROR_READING_VEP_TABLES;
import static no.uib.pap.model.Error.sendError;
import static org.neo4j.driver.v1.Values.parameters;

public class Extractor {

    // Separate sets of objects are created to have more than one attribute of the
    // objects like, iReactions or iPathways
    // Other objects like genes, proteins don't need a separate list, because the
    // identifier is the only attribute used.

    private static ImmutableMap<String, String> iReactions; // Reaction stId to Reaction displayName
    private static ImmutableMap<String, Pathway> iPathways; // Pathway stId to Pathway instance
    private static ImmutableMap<String, String> iProteins; // Protein accession (UniProt) to name
    private static ImmutableSetMultimap<String, String> imapGenesToProteins = null;
    private static ImmutableSetMultimap<String, String> imapEnsemblToProteins = null;
    private static ImmutableSetMultimap<String, String> imapPhysicalEntitiesToReactions = null;
    private static ImmutableSetMultimap<String, String> imapProteinsToReactions = null;
    private static ImmutableSetMultimap<String, String> imapReactionsToPathways = null;
    private static ImmutableSetMultimap<String, String> imapPathwaysToTopLevelPathways = null;
    private static ImmutableSetMultimap<String, Proteoform> imapProteinsToProteoforms = null;
    private static ImmutableSetMultimap<Proteoform, String> imapProteoformsToReactions = null;
    private static ImmutableSetMultimap<String, String> imapRsIdsToProteins = null;
    private static ImmutableSetMultimap<String, String> imapChrBpToProteins = null;
    private static ImmutableSetMultimap<String, String> imapReactionsToParticipants = null;
    private static ImmutableSetMultimap<String, String> imapProteinsToComplexes = null;
    private static ImmutableSetMultimap<String, String> imapComplexesToParticipants = null;

    private static final int rsidColumnIndex = 2;
    // private static final int ensemblColumnIndex = 4;
    private static final int swissprotColumnIndex = 5;
    // private static final int nearestGeneColumnIndex = 7;

    static String outputPath = "../../PathwayMatcher/resources/";

    public static void main(String[] args) {

        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        imapProteinsToReactions = getMapProteinsToReactions();
        System.out.println("Finished map proteins to iReactions.");

        imapRsIdsToProteins = getMapsRsIdsToProteins();
        System.out.println("Finished map rsids to proteins.");

        imapChrBpToProteins = getMapsChrBpToProteins();
        System.out.println("Finished map chrBp to proteins.");

        imapGenesToProteins = getMapGenesToProteins();
        System.out.println("Finished map genes to proteins.");

        imapEnsemblToProteins = getMapEnsemblToProteins();
        System.out.println("Finished map ensembl to proteins.");

        imapProteoformsToReactions = getMapProteoformsToReactions();
        System.out.println("Finished map proteoforms to iReactions.");

        imapReactionsToPathways = getMapReactonsToPathways();
        System.out.println("Finished map iReactions to iPathways.");

        imapPathwaysToTopLevelPathways = getMapPathwaysToTopLevelPathways();
        System.out.println("Finished map iPathways to top level iPathways.");

        iProteins = getProteinNames();
        System.out.println("Finished getting the protein names.");

        imapReactionsToParticipants = getReactionsToParticipants();
        System.out.println("Finished getting reaction participants.");

        imapProteinsToComplexes = getProteinsToComplexes();
        System.out.println("Finished map proteins to complexes.");

        imapComplexesToParticipants = getComplexesToParticipants();
        System.out.println("Finished map of complexes to participants.");

    }

    public static ImmutableSetMultimap<String, String> getComplexesToParticipants() {
        ImmutableSetMultimap.Builder<String, String> builderComplexesToParticipants = new ImmutableSetMultimap.Builder<>();

        //Query the database to fill the data structure
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_COMPLEX_PARTICIPANTS);
        for (Record record : resultList) {
            builderComplexesToParticipants.put(record.get("complex").asString(), record.get("protein").asString());
        }

        imapComplexesToParticipants = builderComplexesToParticipants.build();

        storeSerialized(imapComplexesToParticipants, outputPath + "imapComplexesToParticipants.gz");

        return imapComplexesToParticipants;
    }

    public static ImmutableSetMultimap<String, String> getProteinsToComplexes() {
        ImmutableSetMultimap.Builder<String, String> builderProteinsToComplexes = new ImmutableSetMultimap.Builder<>();

        //Query the database to fill the data structure
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_COMPLEX_PARTICIPANTS);
        for (Record record : resultList) {
            builderProteinsToComplexes.put(record.get("protein").asString(), record.get("complex").asString());
        }

        imapProteinsToComplexes = builderProteinsToComplexes.build();
        storeSerialized(imapProteinsToComplexes, outputPath + "imapProteinsToComplexes.gz");

        return imapProteinsToComplexes;
    }

    public static ImmutableSetMultimap<String, String> getReactionsToParticipants() {
        ImmutableSetMultimap.Builder<String, String> builderReactionsToParticipants = new ImmutableSetMultimap.Builder<>();

        // Query the database and fill the data structure
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_REACTION_PARTICIPANTS);
        for (Record record : resultList) {
            builderReactionsToParticipants.put(record.get("reaction").asString(), record.get("protein").asString());
        }

        imapReactionsToParticipants = builderReactionsToParticipants.build();

        // Serialize list of iReactions
        storeSerialized(imapReactionsToParticipants, outputPath + "imapReactionsToParticipants.gz");

        return imapReactionsToParticipants;
    }

    private static ImmutableSetMultimap<String, String> getMapsRsIdsToProteins() {

        if (imapProteinsToReactions == null) {
            getMapProteinsToReactions();
        }

        ImmutableSetMultimap.Builder<String, String> builderRsIdsToProteins = ImmutableSetMultimap
                .<String, String>builder();

        // Traverse all the vepTables
        for (int chr = 1; chr <= 22; chr++) {
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
        }

        imapRsIdsToProteins = builderRsIdsToProteins.build();

        storeSerialized(imapRsIdsToProteins, outputPath + "imapRsIdsToProteins.gz");

        return imapRsIdsToProteins;
    }

    private static ImmutableSetMultimap<String, String> getMapsChrBpToProteins() {

        if (imapProteinsToReactions == null) {
            getMapProteinsToReactions();
        }

        ImmutableSetMultimap.Builder<String, String> builderChrBpToProteins = ImmutableSetMultimap
                .<String, String>builder();

        // Traverse all the vepTables
        for (int chr = 1; chr <= 22; chr++) {
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
                                builderChrBpToProteins.put(snp.getChr() + "_" + snp.getBp(), protein);
                            }
                        }
                    }
                }
            } catch (IOException ex) {
                sendError(ERROR_READING_VEP_TABLES, chr);
            }
        }

        imapChrBpToProteins = builderChrBpToProteins.build();

        storeSerialized(imapChrBpToProteins, outputPath + "imapChrBpToProteins.gz");
        return imapChrBpToProteins;
    }

    /**
     * Get list of iReactions
     */
    private static ImmutableMap<String, String> getReactions() {

        ImmutableMap.Builder<String, String> builderReactions = ImmutableMap.<String, String>builder();

        // Query the database and fill the data structure
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_ALL_REACTIONS);
        for (Record record : resultList) {
            builderReactions.put(record.get("stId").asString(), record.get("displayName").asString());
        }

        iReactions = builderReactions.build();

        // Serialize list of iReactions
        storeSerialized(iReactions, outputPath + "iReactions.gz");
        return iReactions;
    }

    /**
     * Get list of iPathways
     */
    private static ImmutableMap<String, Pathway> getPathways() {

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

    private static ImmutableSetMultimap<String, String> getMapGenesToProteins() {

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

    private static ImmutableSetMultimap<String, String> getMapEnsemblToProteins() {

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

    private static ImmutableSetMultimap<Proteoform, String> getMapProteoformsToReactions() {
        if (iReactions == null) {
            getReactions();
        }

        // Get mapping from Physical Entities to Reactions
        ImmutableSetMultimap.Builder<String, Proteoform> builderProteinsToProteoforms = ImmutableSetMultimap
                .<String, Proteoform>builder();
        ImmutableSetMultimap.Builder<Proteoform, String> builderProteoformsToReactions = ImmutableSetMultimap
                .<Proteoform, String>builder();
        ImmutableSetMultimap.Builder<String, String> builderPhysicalEntitiesToReactions = ImmutableSetMultimap
                .<String, String>builder();

        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PHYSICALENTITIES_TO_REACTIONS);
        for (Record record : resultList) {
            builderPhysicalEntitiesToReactions.put(record.get("physicalEntity").asString(),
                    record.get("reaction").asString());
        }

        imapPhysicalEntitiesToReactions = builderPhysicalEntitiesToReactions.build();

        // Get mapping from Proteoforms to Physical Entities

        resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PROTEOFORMS_TO_PHYSICALENTITIES);
        for (Record record : resultList) {

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
                for (String reaction : imapPhysicalEntitiesToReactions.get(physicalEntity.toString())) {
                    builderProteoformsToReactions.put(proteoform, reaction);
                }
            }
            builderProteinsToProteoforms.put(record.get("protein").asString(), proteoform);
        }

        imapProteinsToProteoforms = builderProteinsToProteoforms.build();
        imapProteoformsToReactions = builderProteoformsToReactions.build();

        storeSerialized(imapProteinsToProteoforms, outputPath + "imapProteinsToProteoforms.gz");
        storeSerialized(imapProteoformsToReactions, outputPath + "imapProteoformsToReactions.gz");

        return imapProteoformsToReactions;
    }

    /**
     * Get mapping from proteins to iReactions
     */
    private static ImmutableSetMultimap<String, String> getMapProteinsToReactions() {

        if (iReactions == null) {
            getReactions();
        }

        // Query the database and fill the data structure
        ImmutableSetMultimap.Builder<String, String> builderProteinsToReactions = ImmutableSetMultimap
                .<String, String>builder();
        List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PROTEINS_TO_REACTIONS);

        for (Record record : resultList) {
            builderProteinsToReactions.put(record.get("protein").asString(), record.get("reaction").asString());
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
     * Get mapping from iReactions to iPathways
     */
    private static ImmutableSetMultimap<String, String> getMapReactonsToPathways() {

        if (iReactions == null) {
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
    private static ImmutableSetMultimap<String, String> getMapPathwaysToTopLevelPathways() {

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
        }
    }
}
