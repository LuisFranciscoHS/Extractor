package no.uib.pap.extractor.neo4j;

import static no.uib.pap.model.Error.ERROR_READING_VEP_TABLES;
import static no.uib.pap.model.Error.sendError;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.io.Reader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang3.tuple.Pair;
import org.neo4j.driver.v1.Record;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.LinkedListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeBasedTable;
import com.google.common.collect.TreeMultimap;
import com.google.common.collect.ImmutableSetMultimap.Builder;
import com.google.common.collect.ImmutableTable;

import no.uib.pap.model.Pathway;
import no.uib.pap.model.Proteoform;
import no.uib.pap.model.ProteoformFormat;
import no.uib.pap.model.Snp;

public class Extractor {

	// Separate sets of objects are created to have more than one attribute of the
	// objects like, reactions or pathways
	// Other objects like genes, proteins don't need a separate list, because the
	// identifier is the only attribute used.

	private static ImmutableMap<String, String> reactions; // Reaction stId to Reaction displayName
	private static ImmutableMap<String, Pathway> pathways; // Pathway stId to Pathway instance
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

	private static final int rsidColumnIndex = 2;
	// private static final int ensemblColumnIndex = 4;
	private static final int swissprotColumnIndex = 5;
	// private static final int nearestGeneColumnIndex = 7;

	public static void main(String[] args) {

		ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

		getMapProteinsToReactions();
		System.out.println("Finished map proteins to reactions.");

		getMapsSnpsToProteins();
		System.out.println("Finished maps snps to proteins.");

		getMapGenesToProteins();
		System.out.println("Finished map genes to proteins.");

		getMapEnsemblToProteins();
		System.out.println("Finished map ensembl to proteins.");

		getMapProteoformsToReactions();
		System.out.println("Finished map proteoforms to reactions.");

		getMapReactonsToPathways();
		System.out.println("Finished map reactions to pathways.");

		getMapPathwaysToTopLevelPathways();
		System.out.println("Finished map pathways to top level pathways.");
	}

	private static void getMapsSnpsToProteins() {

		if (imapProteinsToReactions == null) {
			getMapProteinsToReactions();
		}

		ImmutableSetMultimap.Builder<String, String> builderRsIdsToProteins = ImmutableSetMultimap
				.<String, String>builder();

		ImmutableSetMultimap.Builder<String, String> builderChrBpToProteins = ImmutableSetMultimap
				.<String, String>builder();

		// Traverse all the vepTables
		for (int chr = 1; chr <= 22; chr++) {
			System.out.println("Scanning vepTable for chromosome " + chr);
			try {
				BufferedReader br = getBufferedReaderFromResource(chr + ".gz");
				br.readLine(); // Read header line

				for (String line; (line = br.readLine()) != null;) {

					Multimap<Snp, String> snpToSwissprotMap = getSNPAndSwissProtFromVep(line);

					for (Map.Entry<Snp, String> snpToSwissprotPair : snpToSwissprotMap.entries()) {
						Snp snp = snpToSwissprotPair.getKey();
						String protein = snpToSwissprotPair.getValue();

						if (!protein.equals("NA")) {
							if (imapProteinsToReactions.containsKey(protein)) {
								builderRsIdsToProteins.put(snp.getRsid(), protein);
								builderChrBpToProteins.put(snp.getChr() + "_" + snp.getBp(), protein);
							}
						}
					}
				}
			} catch (IOException ex) {
				sendError(ERROR_READING_VEP_TABLES, chr);
			}
		}

		imapRsIdsToProteins = builderRsIdsToProteins.build();
		imapChrBpToProteins = builderChrBpToProteins.build();

		storeSerialized(imapChrBpToProteins, "imapChrBpToProteins.gz");
		storeSerialized(imapRsIdsToProteins, "imapRsIdsToProteins.gz");
	}

	/**
	 * Get list of reactions
	 */
	private static void getReactions() {

		ImmutableMap.Builder<String, String> builderReactions = ImmutableMap.<String, String>builder();

		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_ALL_REACTIONS);
		for (Record record : resultList) {
			builderReactions.put(record.get("stId").asString(), record.get("displayName").asString());
		}

		reactions = builderReactions.build();

		// Serialize list of reactions
		storeSerialized(reactions, "reactions.gz");
	}

	/**
	 * Get list of pathways
	 */
	private static void getPathways() {

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

		pathways = builderPathways.build();

		// Serialize pathways
		storeSerialized(pathways, "pathways.gz");
	}

	private static void getMapGenesToProteins() {

		// Query the database and fill the data structure
		ImmutableSetMultimap.Builder<String, String> builderGenesToProteins = ImmutableSetMultimap
				.<String, String>builder();
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_GENES_TO_PROTEINS);
		for (Record record : resultList) {
			builderGenesToProteins.put(record.get("gene").asString(), record.get("protein").asString());
		}
		imapGenesToProteins = builderGenesToProteins.build();

		storeSerialized(imapGenesToProteins, "imapGenesToProteins.gz");
	}

	private static void getMapEnsemblToProteins() {

		// Query the database and fill the data structure
		ImmutableSetMultimap.Builder<String, String> builderEnsemblToProteins = ImmutableSetMultimap
				.<String, String>builder();
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_ENSEMBL_TO_PROTEINS);
		for (Record record : resultList) {
			builderEnsemblToProteins.put(record.get("ensembl").asString(), record.get("uniprot").asString());
		}
		imapEnsemblToProteins = builderEnsemblToProteins.build();

		storeSerialized(imapEnsemblToProteins, "imapEnsemblToProteins.gz");
	}

	private static void getMapProteoformsToReactions() {
		if (reactions == null) {
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
			TreeMultimap<String, Long> ptms = TreeMultimap.create();

			for (Object ptm : record.get("ptms").asList()) {
				String[] parts = ptm.toString().split(":");
				String mod = parts[0];
				Long site = Proteoform.interpretCoordinateFromStringToLong(parts[1]);
				ptms.put(mod, site);
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

		storeSerialized(imapProteinsToProteoforms, "imapProteinsToProteoforms.gz");
		storeSerialized(imapProteoformsToReactions, "imapProteoformsToReactions.gz");
	}

	/**
	 * Get mapping from proteins to reactions
	 */
	private static void getMapProteinsToReactions() {

		if (reactions == null) {
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

		storeSerialized(imapProteinsToReactions, "imapProteinsToReactions.gz");
	}

	/**
	 * Get mapping from reactions to pathways
	 */
	private static void getMapReactonsToPathways() {

		if (reactions == null) {
			getReactions();
		}
		if (pathways == null) {
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

		storeSerialized(imapReactionsToPathways, "imapReactionsToPathways.gz");
	}

	/**
	 * Get mapping from pathways to top level pathways
	 */
	private static void getMapPathwaysToTopLevelPathways() {

		if (pathways.size() == 0) {
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

		storeSerialized(imapPathwaysToTopLevelPathways, "imapPathwaysToTopLevelPathways.gz");
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
