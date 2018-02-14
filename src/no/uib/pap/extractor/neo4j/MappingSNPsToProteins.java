package no.uib.pap.extractor.neo4j;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.HashSet;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import org.neo4j.driver.v1.Record;
import org.neo4j.driver.v1.Session;
import org.neo4j.driver.v1.StatementResult;

import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;

import no.uib.pap.model.Snp;

public class MappingSNPsToProteins {
	private static final int rsidVepColumn = 2;
	private static final int swissprotVepColumn = 5;

	static void getMappingSNPsToProteins() {
		HashSet<String> proteinSet = new HashSet<>(16000); // Lowest power of 2 that can hold all the proteins ~10500
		TreeMultimap<Snp, String> allSnpToSwissprotMap = TreeMultimap.create();

		System.out.println("Extracting map from SNPs to Proteins...");

		// Get list of proteins available in Reactome
		try {
			Session session = ConnectionNeo4j.driver.session();

			String query = ReactomeQueries.GET_ALL_PROTEINS;
			StatementResult queryResult;

			queryResult = session.run(query);

			while (queryResult.hasNext()) {
				Record record = queryResult.next();
				String uniprotAccession = record.get("protein").asString();
				proteinSet.add(uniprotAccession);
			}

			session.close();
		} catch (org.neo4j.driver.v1.exceptions.ClientException ex) {
			System.err.println("ClientException: " + ex.getMessage());
		}

		
		
		// Traverse VEP tables of each chromosome
		for (int chr = 21; chr <= 22; chr++) {
			System.out.println("Scanning vepTable for chromosome " + chr);
			try {
				BufferedReader br = no.uib.pap.utils.FileUtils.getBufferedReader("resources/VEP/" + chr + ".gz");
				br.readLine(); // Read header line

				for (String line; (line = br.readLine()) != null;) {

					Multimap<Snp, String> snpToSwissprotMap = getSNPAndSwissProtFromVep(line);

					for (Map.Entry<Snp, String> snpToSwissprotPair : snpToSwissprotMap.entries()) {
						if (!snpToSwissprotPair.getValue().equals("NA")) {
							if(proteinSet.contains(snpToSwissprotPair.getValue())) {
								allSnpToSwissprotMap.put(snpToSwissprotPair.getKey(), snpToSwissprotPair.getValue());
							}
						}
						else {
							break;
						}
					}
				}
			} catch (IOException ex) {
				System.err.println("IOException: " + ex.getMessage());
			}
		}
		
		System.out.println("Serializing map...");

		FileOutputStream fos = null;
		GZIPOutputStream gz;
		ObjectOutputStream oos;
		try {
			fos = new FileOutputStream("snpsToProteins.gz");
			gz = new GZIPOutputStream(fos);
			oos = new ObjectOutputStream(gz);
			oos.writeObject(allSnpToSwissprotMap);
			oos.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.out.println("Finished map from SNPs to Proteins");
	}

	public static Multimap<Snp, String> getSNPAndSwissProtFromVep(String line) {
		TreeMultimap<Snp, String> mapping = TreeMultimap.create();
		String[] fields = line.split(" ");
		Integer chr = Integer.valueOf(fields[0]);
		Long bp = Long.valueOf(fields[1]);

		String[] rsids = fields[rsidVepColumn].split(",");
		String[] uniprots = fields[swissprotVepColumn].split(",");

		for (String rsid : rsids) {
			for (String uniprot : uniprots) {
				Snp snp = new Snp(chr, bp, rsid);
				mapping.put(snp, uniprot);
			}
		}
		return mapping;
	}
}
