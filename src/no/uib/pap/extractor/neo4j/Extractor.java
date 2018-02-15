package no.uib.pap.extractor.neo4j;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.neo4j.driver.v1.Record;

import com.google.common.collect.LinkedListMultimap;
import com.google.common.collect.TreeMultimap;

import no.uib.pap.model.Pathway;
import no.uib.pap.model.Proteoform;
import no.uib.pap.model.ProteoformFormat;

public class Extractor {

	// Separate sets of objects are created to have more than one attribute of the
	// objects like, reactions or pathways
	// Other objects like genes, proteins don't need a separate list, because the
	// identifier is the only attribute used.

	private static HashMap<String, String> reactions = new HashMap<>(11500); // Reaction stId to Reaction displayName
	private static HashMap<String, Pathway> pathways = new HashMap<>(); // Pathway stId to Pathway instance
	private static TreeMultimap<String, String> mapGenesToProteins = TreeMultimap.create(); // Use multimap because
	private static TreeMultimap<String, String> mapEnsemblToProteins = TreeMultimap.create();
	private static TreeMultimap<String, String> mapPhysicalEntitiesToReactions = TreeMultimap.create();
	private static TreeMultimap<String, String> mapProteinsToReactions = TreeMultimap.create();
	private static TreeMultimap<String, String> mapReactionsToPathways = TreeMultimap.create();
	private static TreeMultimap<String, String> mapPathwaysToTopLevelPathways = TreeMultimap.create();
	private static TreeMultimap<String, Proteoform> mapProteinsToProteoforms = TreeMultimap.create();
	private static TreeMultimap<Proteoform, String> mapProteoformsToReactions = TreeMultimap.create();
	// TODO Change these structures immutable structures. 

	public static void main(String[] args) {

		ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

		getMapGenesToProteins();
		System.out.println("Finished map genes to proteins.");
		
		getMapEnsemblToProteins();
		System.out.println("Finished map ensembl to proteins.");
		
		getMapProteinsToReactions();
		System.out.println("Finished map proteins to reactions.");

		getMapProteoformsToReactions();
		System.out.println("Finished map proteoforms to reactions.");

		getMapReactonsToPathways();
		System.out.println("Finished map reactions to pathways.");

		getMapPathwaysToTopLevelPathways();
		System.out.println("Finished map pathways to top level pathways.");

	}

	/**
	 * Get list of reactions
	 */
	private static void getReactions() {

		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_ALL_REACTIONS);
		for (Record record : resultList) {
			reactions.put(record.get("stId").asString(), record.get("displayName").asString());
		}

		// Serialize list of reactions
		storeSerialized(reactions, "reactions.gz");
	}

	/**
	 * Get list of pathways
	 */
	private static void getPathways() {

		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_ALL_PATHWAYS);
		for (Record record : resultList) {
			
			// Get first the attributes that can be obtained directly
			Pathway pathway = new Pathway(record.get("stId").asString(), record.get("displayName").asString());
			pathway.setNumEntitiesTotal(record.get("numEntitiesTotal").asInt());
			pathway.setNumReactionsTotal(record.get("numReactionsTotal").asInt());
			pathway.setNumReactionsTotal(record.get("numProteoformsTotal").asInt());
			pathways.put(pathway.getStId(), pathway);
		}

		// Serialize pathways
		storeSerialized(pathways, "pathways.gz");
	}

	private static void getMapGenesToProteins() {
		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_GENES_TO_PROTEINS);
		for (Record record : resultList) {
			mapGenesToProteins.put(record.get("gene").asString(), record.get("protein").asString());
		}
		
		storeSerialized(mapGenesToProteins, "mapGenesToProteins.gz");
	}

	private static void getMapEnsemblToProteins() {
		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_ENSEMBL_TO_PROTEINS);
		for (Record record : resultList) {
			mapEnsemblToProteins.put(record.get("ensembl").asString(), record.get("uniprot").asString());
		}
		
		storeSerialized(mapEnsemblToProteins, "mapEnsemblToProteins.gz");
	}

	private static void getMapProteoformsToReactions() {
		if (reactions.size() == 0) {
			getReactions();
		}

		// Get mapping from Physical Entities to Reactions
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PHYSICALENTITIES_TO_REACTIONS);
		for (Record record : resultList) {
			mapPhysicalEntitiesToReactions.put(
					record.get("physicalEntity").asString(),
					record.get("reaction").asString());
		}

		// Get mapping from Proteoforms to Physical Entities

		resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PROTEOFORMS_TO_PHYSICALENTITIES);
		for (Record record : resultList) {
						
			String isoform = record.get("isoform").asString();
			TreeMultimap<String, Long> ptms = TreeMultimap.create();

			for(Object ptm : record.get("ptms").asList()) {
				String[] parts = ptm.toString().split(":");
				String mod = parts[0];
				Long site = Proteoform.interpretCoordinateFromStringToLong(parts[1]);
				ptms.put(mod, site);
			}
			
			Proteoform proteoform = new Proteoform(isoform, ptms);

			for(Object physicalEntity : record.get("peSet").asList()) {
				for(String reaction : mapPhysicalEntitiesToReactions.get(physicalEntity.toString())) {
					mapProteoformsToReactions.put(proteoform, reaction);
				}
			}
			
			mapProteinsToProteoforms.put(record.get("protein").asString(), proteoform);
		}

		storeSerialized(mapProteinsToProteoforms, "mapProteinsToProteoforms.gz");
		storeSerialized(mapProteoformsToReactions, "mapProteoformsToReactions.gz");
	}

	/**
	 * Get mapping from proteins to reactions
	 */
	private static void getMapProteinsToReactions() {

		if (reactions.size() == 0) {
			getReactions();
		}

		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_PROTEINS_TO_REACTIONS);
		for (Record record : resultList) {
			mapProteinsToReactions.put(record.get("protein").asString(), record.get("reaction").asString());
		}

		storeSerialized(mapProteinsToReactions, "mapProteinsToReactions.gz");
	}

	/**
	 * Get mapping from reactions to pathways
	 */
	private static void getMapReactonsToPathways() {

		if (reactions.size() == 0) {
			getReactions();
		}
		if (pathways.size() == 0) {
			getPathways();
		}

		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GET_MAP_REACTIONS_TO_PATHWAYS);
		for (Record record : resultList) {
			mapReactionsToPathways.put(record.get("reaction").asString(), record.get("pathway").asString());
		}

		storeSerialized(mapReactionsToPathways, "mapReactionsToPathways.gz");
	}

	/**
	 * Get mapping from pathways to top level pathways
	 */
	private static void getMapPathwaysToTopLevelPathways() {

		if (pathways.size() == 0) {
			getPathways();
		}

		// Query the database and fill the data structure
		List<Record> resultList = ConnectionNeo4j.query(ReactomeQueries.GEP_MAP_PATHWAYS_TO_TOPLEVELPATHWAYS);
		for (Record record : resultList) {
			mapReactionsToPathways.put(record.get("reaction").asString(), record.get("pathway").asString());
		}

		storeSerialized(mapPathwaysToTopLevelPathways, "mapPathwaysToTopLevelPathways.gz");

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
