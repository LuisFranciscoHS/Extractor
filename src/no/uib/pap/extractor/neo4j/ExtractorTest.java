package no.uib.pap.extractor.neo4j;

import com.google.common.collect.ImmutableMap;
import no.uib.pap.model.Pathway;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.v1.Record;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class ExtractorTest {

    @org.junit.jupiter.api.BeforeEach
    void setUp() {

    }

    @Test
    void genesToProteinsTest() {
//
//        for(String geneName : imapGenesToProteins.keySet()) {
//            System.out.println(geneName);
//        }
//
//        String geneName = "INSR";
//        System.out.print(geneName + ": ");
//        for(String protein : imapGenesToProteins.get(geneName)) {
//            System.out.print(protein + ", ");
//        }
//        System.out.println("");
    }

    @Test
    void test() {
//         for(Pathway pathway : iPathways){
//         System.out.println(pathway.toString());
//         }
//
//         for (String reaction : imapProteinsToReactions.get("P01308")) {
//         for (String pathway : imapReactionsToPathways.get(reaction)) {
//         System.out.println(pathway + "\t" + iPathways.get(pathway).getDisplayName() +
//         "\t" + reaction + "\t"
//         + iReactions.get(reaction));
//         }
//         }
    }

    @Test
    void getProteinNamesTest(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableMap<String,String> iProteins = Extractor.getProteinNames();

        assertEquals("Hemoglobin subunit beta", iProteins.get("P68871"));
        assertEquals("Hemoglobin subunit alpha", iProteins.get("P69905"));
        assertEquals("Insulin", iProteins.get("P01308"));
    }


    @Test
    void getConnectionsMap() {

    }

    @Test
    void getSNPAndSwissProtFromVep() {
    }
}