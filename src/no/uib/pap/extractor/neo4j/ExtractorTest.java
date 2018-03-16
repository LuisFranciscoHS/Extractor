package no.uib.pap.extractor.neo4j;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSetMultimap;
import no.uib.pap.model.Reaction;
import no.uib.pap.model.Role;
import org.junit.jupiter.api.Test;

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
    void getEnsembleToProteinsTest(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");
        ImmutableSetMultimap<String,String> imapEnsembleToProteins = Extractor.getMapEnsemblToProteins();

        assertEquals(1, imapEnsembleToProteins.get("ENSP00000380389").size());
        assertTrue(imapEnsembleToProteins.get("ENSP00000380389").contains("P19438"));

        assertEquals(1, imapEnsembleToProteins.get("ENSG00000143226").size());
        assertTrue(imapEnsembleToProteins.get("ENSG00000143226").contains("P12318"));
    }

    @Test
    void getProteinNamesTest() {
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableMap<String, String> iProteins = Extractor.getProteinNames();

        assertEquals("Hemoglobin subunit beta", iProteins.get("P68871"));
        assertEquals("Hemoglobin subunit alpha", iProteins.get("P69905"));
        assertEquals("Insulin", iProteins.get("P01308"));
    }

    @Test
    void getReactionsNeighboursTest1(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");
        ImmutableMap<String,Reaction> iReactions = Extractor.getReactions();

        assertEquals(4, iReactions.get("R-HSA-2230938").getParticipants().keySet().size());
        assertTrue(iReactions.get("R-HSA-2230938").getParticipants().containsKey("P00738"));
        assertTrue(iReactions.get("R-HSA-2230938").getParticipants().containsKey("P69905"));
        assertTrue(iReactions.get("R-HSA-2230938").getParticipants().containsKey("Q86VB7"));
        assertTrue(iReactions.get("R-HSA-2230938").getParticipants().containsKey("P68871"));
    }

    @Test
    void getReactionNeightboursTest2() {
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableMap<String, Reaction> iReactions = Extractor.getReactions();

        assertEquals(2, iReactions.get("R-HSA-74716").getParticipants().keySet().size());
        assertTrue(iReactions.get("R-HSA-74716").getParticipants().containsKey("P06213"));
        assertTrue(iReactions.get("R-HSA-74716").getParticipants().containsKey("P01308"));

        assertEquals(1, iReactions.get("R-HSA-74730").getParticipants().keySet().size());
        assertTrue(iReactions.get("R-HSA-74730").getParticipants().containsKey("P01308"));

        assertEquals(2, iReactions.get("R-HSA-74716").getParticipants().keySet().size());
        assertTrue(iReactions.get("R-HSA-74716").getParticipants().containsKey("P01308"));
        assertTrue(iReactions.get("R-HSA-74716").getParticipants().containsKey("P06213"));

        assertEquals(2, iReactions.get("R-HSA-74726").getParticipants().keySet().size());
        assertTrue(iReactions.get("R-HSA-74726").getParticipants().containsKey("P01308"));
        assertTrue(iReactions.get("R-HSA-74726").getParticipants().containsKey("P06213"));
    }

    @Test
    void getReactionParticipantRolesTest(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableMap<String, Reaction> iReactions = Extractor.getReactions();

        assertEquals(14, iReactions.get("R-HSA-8863895").getParticipants().keySet().size());
        assertTrue(iReactions.get("R-HSA-8863895").getParticipants().get("O14920").contains(Role.CATALYSTACTIVITY));
        assertTrue(iReactions.get("R-HSA-8863895").getParticipants().get("O00161").contains(Role.INPUT));
        assertTrue(iReactions.get("R-HSA-8863895").getParticipants().get("O00161").contains(Role.OUTPUT));
        assertTrue(iReactions.get("R-HSA-8863895").getParticipants().get("Q99836").contains(Role.REGULATEDBY));

    }

    @Test
    void imapProteinsToComplexesTest() {
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapProteinsToComplexes = Extractor.getProteinsToComplexes();

        assertEquals(32, imapProteinsToComplexes.get("Q9Y297").size());
        assertTrue(imapProteinsToComplexes.get("Q9Y297").contains("R-HSA-174138"));
        assertTrue(imapProteinsToComplexes.get("Q9Y297").contains("R-HSA-8952593"));

        assertEquals(1, imapProteinsToComplexes.get("P00558").size());
        assertTrue(imapProteinsToComplexes.get("P00558").contains("R-HSA-70484"));
    }

    @Test
    void imapComplexesToComponentsTest() {
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapComplexesToParticipants = Extractor.getComplexesToComponents();

        assertEquals(1, imapComplexesToParticipants.get("R-HSA-70484").size());
        assertTrue(imapComplexesToParticipants.get("R-HSA-70484").contains("P00558"));

        assertEquals(6, imapComplexesToParticipants.get("R-HSA-174138").size());
        assertTrue(imapComplexesToParticipants.get("R-HSA-174138").contains("Q9Y297"));
        assertTrue(imapComplexesToParticipants.get("R-HSA-174138").contains("P63208"));
        assertTrue(imapComplexesToParticipants.get("R-HSA-174138").contains("Q9UKT4"));
    }

    @Test
    void imapSetsToMembersAndCandidatesTest() {
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapSetsToMembersAndCandidates = Extractor.getSetsToMembersAndCandidates();

        assertEquals(3, imapSetsToMembersAndCandidates.get("R-HSA-1008234").size());
        assertTrue(imapSetsToMembersAndCandidates.get("R-HSA-1008234").contains("Q9ULX9"));

        assertEquals(3, imapSetsToMembersAndCandidates.get("R-HSA-114528").size());
        assertTrue(imapSetsToMembersAndCandidates.get("R-HSA-114528").contains("P17252"));
        assertTrue(imapSetsToMembersAndCandidates.get("R-HSA-114528").contains("P05771"));
        assertTrue(imapSetsToMembersAndCandidates.get("R-HSA-114528").contains("P05129"));
    }

    @Test
    void imapRsidsToProteinsTest(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapRsidsToProteins = Extractor.getRsIdsToProteins(11);

        assertTrue(imapRsidsToProteins.containsKey("rs10840447"));
        assertTrue(imapRsidsToProteins.containsKey("rs7110099"));
        assertTrue(imapRsidsToProteins.containsKey("rs555583938"));
        assertTrue(imapRsidsToProteins.get("rs10840447").contains("P01308"));
        assertTrue(imapRsidsToProteins.get("rs7110099").contains("P01308"));
        assertTrue(imapRsidsToProteins.get("rs555583938").contains("P01308"));
    }

    @Test
    void imapChrBpToProteinsTest(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<Long,String> imapChrBpToProteins = Extractor.getChrBpToProteins(11);

        assertTrue(imapChrBpToProteins.containsKey(2176042L));
        assertTrue(imapChrBpToProteins.containsKey(2176105L));
        assertTrue(imapChrBpToProteins.containsKey(2176134L));
        assertTrue(imapChrBpToProteins.get(2176042L).contains("P01308"));
        assertTrue(imapChrBpToProteins.get(2176105L).contains("P01308"));
        assertTrue(imapChrBpToProteins.get(2176134L).contains("P01308"));
    }

    @Test
    void imapGeneticVariantsToProteinsTest5(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapRsidsToProteins = Extractor.getRsIdsToProteins(5);

        assertTrue(imapRsidsToProteins.containsKey("rs17238540"));
        assertTrue(imapRsidsToProteins.get("rs17238540").contains("P04035"));

        ImmutableSetMultimap<Long,String> imapChrBpToProteins = Extractor.getChrBpToProteins(5);

        assertTrue(imapChrBpToProteins.containsKey(74655498L));
        assertTrue(imapChrBpToProteins.get(74655498L).contains("P04035"));
    }

    @Test
    void imapGeneticVariantsToProteinsTest19(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapRsidsToProteins = Extractor.getRsIdsToProteins(19);

        assertTrue(imapRsidsToProteins.containsKey("rs12979860"));
        assertTrue(imapRsidsToProteins.containsKey("rs12980275"));
        assertTrue(imapRsidsToProteins.get("rs12979860").contains("Q8IZI9"));
        assertTrue(imapRsidsToProteins.get("rs12979860").contains("Q8IZI9"));

        ImmutableSetMultimap<Long,String> imapChrBpToProteins = Extractor.getChrBpToProteins(19);

        assertTrue(imapChrBpToProteins.containsKey(39729266L));
        assertTrue(imapChrBpToProteins.containsKey(39729326L));
        assertTrue(imapChrBpToProteins.get(39729266L).contains("Q8IZI9"));
        assertTrue(imapChrBpToProteins.get(39729326L).contains("Q8IZI9"));
    }

    @Test
    void imapGeneticVariantsToProteinsTest1(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapRsidsToProteins = Extractor.getRsIdsToProteins(1);

        assertTrue(imapRsidsToProteins.containsKey("rs2816958"));
        assertTrue(imapRsidsToProteins.get("rs2816958").contains("O00482"));

        ImmutableSetMultimap<Long,String> imapChrBpToProteins = Extractor.getChrBpToProteins(1);

        assertTrue(imapChrBpToProteins.containsKey(200101920L));
        assertTrue(imapChrBpToProteins.get(200101920L).contains("O00482"));
    }

    @Test
    void imapGeneticVariantsToProteinsTest22(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapRsidsToProteins = Extractor.getRsIdsToProteins(22);

        assertTrue(imapRsidsToProteins.containsKey("rs35873774"));
        assertTrue(imapRsidsToProteins.get("rs35873774").contains("P17861"));

        ImmutableSetMultimap<Long,String> imapChrBpToProteins = Extractor.getChrBpToProteins(22);

        assertTrue(imapChrBpToProteins.containsKey(29191932L));
        assertTrue(imapChrBpToProteins.get(29191932L).contains("P17861"));
    }

    @Test
    void imapGeneticVariantsToProteinsTest20(){
        ConnectionNeo4j.initializeNeo4j("bolt://127.0.0.1:7687", "", "");

        ImmutableSetMultimap<String,String> imapRsidsToProteins = Extractor.getRsIdsToProteins(20);

        assertTrue(imapRsidsToProteins.containsKey("rs1883832"));
        assertTrue(imapRsidsToProteins.get("rs1883832").contains("P25942"));

        ImmutableSetMultimap<Long,String> imapChrBpToProteins = Extractor.getChrBpToProteins(20);

        assertTrue(imapChrBpToProteins.containsKey(44746982L));
        assertTrue(imapChrBpToProteins.get(44746982L).contains("P25942"));
    }

    @Test
    void getConnectionsMap() {

    }

    @Test
    void getSNPAndSwissProtFromVep() {
    }
}