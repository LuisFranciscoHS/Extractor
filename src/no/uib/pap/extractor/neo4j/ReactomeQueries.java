package no.uib.pap.extractor.neo4j;

public interface ReactomeQueries {
    public static final String GET_MAP_GENES_TO_PROTEINS = "MATCH (ewas:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nWITH re.identifier as protein, re.geneName as genes\nWHERE size(genes) > 0  \nUNWIND genes as gene\nRETURN DISTINCT gene, protein";

    public static final String GET_MAP_ENSEMBL_TO_PROTEINS = "MATCH (ewas:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nUNWIND re.otherIdentifier as ensembl\nWITH DISTINCT ensembl, re.identifier as uniprot\nWHERE ensembl STARTS WITH \"ENS\"\nRETURN ensembl, uniprot";

    public static final String GET_ALL_PROTEINS = "MATCH (pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\nRETURN DISTINCT re.identifier as protein";

    public static final String GET_COUNT_ALL_PROTEINS = "MATCH (pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\nRETURN count(DISTINCT re.identifier) as count";

    public static final String GET_ALL_PROTEINS_WITH_ISOFORMS = "MATCH (pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\n" +
            "RETURN DISTINCT (CASE WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier ELSE re.identifier END) as Identifiers";

    public static final String GET_COUNT_ALL_PROTEINS_WITH_ISOFORMS = "MATCH (pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\n" +
            "RETURN count(DISTINCT CASE WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier ELSE re.identifier END) as count";

    public static final String GET_ALL_REACTIONS = "MATCH (n:ReactionLikeEvent{speciesName:\"Homo sapiens\"}) \nRETURN DISTINCT n.stId as stId, n.displayName as displayName";

    public static final String GET_COUNT_ALL_REACTIONS = "MATCH (n:ReactionLikeEvent{speciesName:\"Homo sapiens\"}) RETURN count(DISTINCT n) as count";

    public static final String GET_MAP_PROTEINS_TO_REACTIONS = "MATCH (r:ReactionLikeEvent{speciesName:\"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\nRETURN DISTINCT re.identifier as protein, r.stId AS reaction";

    public static final String GET_MAP_PHYSICALENTITIES_TO_REACTIONS = "MATCH (r:ReactionLikeEvent{speciesName:\"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"})\nRETURN DISTINCT pe.stId as physicalEntity, r.stId AS reaction";

    public static final String GET_ALL_PATHWAYS = "MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nWITH DISTINCT  pe, re\nOPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)\nWITH DISTINCT \n  pe,\n  re.identifier AS protein,\n  re.variantIdentifier AS isoform,\n  tm.coordinate as coordinate, \n  mod.identifier as type ORDER BY type, coordinate\nWITH DISTINCT   \n  pe,\n  protein,\n  isoform,\n  COLLECT(CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE \"null\" END + \":\" + type) AS ptms \n WITH DISTINCT   \n  pe,\n  protein,\n  COLLECT(CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END + ptms) as proteoform\n  MATCH (p:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName: \"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe)\nRETURN DISTINCT p.stId AS stId, p.displayName as displayName, count(DISTINCT rle) as numReactionsTotal, count(DISTINCT protein) as numEntitiesTotal, count(DISTINCT proteoform) as numProteoformsTotal\n  ORDER BY stId";

    public static final String GET_COUNT_ALL_PATHWAYS = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"}) RETURN count(DISTINCT p) as count";

    public static final String GET_MAP_REACTIONS_TO_PATHWAYS = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName: \"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\nRETURN DISTINCT p.stId AS pathway, rle.stId as reaction";

    public static final String GEP_MAP_PATHWAYS_TO_TOPLEVELPATHWAYS = "MATCH (tlp:TopLevelPathway{speciesName:'Homo sapiens'})-[:hasEvent*]->(p:Pathway{speciesName:'Homo sapiens'})\n" +
            "RETURN DISTINCT p.stId AS pathway, tlp.stId as topLevelPathway";

    public static final String GET_MAP_PROTEOFORMS_TO_PHYSICALENTITIES = "MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nWITH DISTINCT pe, re\nOPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)\nWITH DISTINCT pe.stId AS physicalEntity,\n                re.identifier AS protein,\n                re.variantIdentifier AS isoform,\n                tm.coordinate as coordinate, \n                mod.identifier as type ORDER BY type, coordinate\nWITH DISTINCT physicalEntity,\n\t\t\t\tprotein,\n                CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END as isoform,\n                COLLECT(type + \":\" + CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE \"null\" END) AS ptms\n                RETURN protein, isoform, ptms, collect(physicalEntity) as peSet\n                ORDER BY isoform, ptms";

    static final String GET_HIT_COUNT_BY_PROTEIN_PROTEOFORMS = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName: \"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt', identifier:\"P31749\"})\nWITH DISTINCT  p, rle, pe, re\nOPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)\nWITH DISTINCT \n  p.stId as pathway,\n  rle.stId as reaction,\n  pe.stId as pe,\n  re.identifier AS protein,\n  re.variantIdentifier AS isoform,\n  tm.coordinate as coordinate, \n  mod.identifier as type ORDER BY type, coordinate\nWITH DISTINCT  pathway, reaction, \n  pe,\n  protein,\n  isoform,\n  COLLECT(CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE \"null\" END + \":\" + type) AS ptms \n WITH DISTINCT pathway, reaction, protein, COLLECT(CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END + ptms) as proteoform\n  RETURN  DISTINCT   count(pathway), count(reaction), protein, proteoform\n  ORDER BY proteoform";

    static final String GET_MAPPING_BY_PROTEIN_LIST = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName:\"Homo sapiens\"}),\n      (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"}),\n      (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\n      WHERE re.identifier IN [\"P01308\"]\nRETURN DISTINCT re.identifier, rle.stId, rle.displayName, p.stId, p.displayName\nORDER BY rle.stId";

    static final String GET_MAPPING_BY_PROTEOFORMS_BY_PROTEIN = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName: \"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\n" +
            "WHERE re.identifier = \"P01308\"\n" +
            "WITH DISTINCT p, rle, pe, re\n" +
            "OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)\n" +
            "WITH DISTINCT p, rle, pe.stId AS physicalEntity,\n" +
            "                re.identifier AS protein,\n" +
            "                re.variantIdentifier AS isoform,\n" +
            "                tm.coordinate as coordinate, \n" +
            "                mod.identifier as type ORDER BY type, coordinate\n" +
            "WITH DISTINCT p, rle, physicalEntity,\n" +
            "\t\t\t\tprotein,\n" +
            "                CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END as isoform,\n" +
            "                COLLECT(type + \":\" + CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE \"null\" END) AS ptms\n" +
            "RETURN DISTINCT p.stId as pathway, rle.stId as reaction, (CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END + ptms) as proteoform       \n" +
            "                ORDER BY proteoform";

    static final String GET_PROTEOFORMS_BY_PATHWAY = "MATCH (pathway:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName:\"Homo sapiens\"}),\n" +
            "      (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\n" +
            "      WHERE pathway.stId IN [\"R-HSA-5683057\",\"R-HSA-162582\"]\n" +
            "WITH DISTINCT pathway, rle, pe, re\n" +
            "OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)\n" +
            "WITH DISTINCT pathway, rle.stId as reaction, pe.stId AS physicalEntity,\n" +
            "                re.identifier AS protein, re.variantIdentifier AS isoform,  tm.coordinate as coordinate, \n" +
            "                mod.identifier as type \n" +
            "ORDER BY type, coordinate\n" +
            "WITH DISTINCT pathway, reaction, physicalEntity, protein,\n" +
            "                CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END as isoform,\n" +
            "                COLLECT(type + \":\" + CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE \"null\" END) AS ptms\n" +
            "RETURN pathway.stId, pathway.displayName, reaction, protein, isoform, ptms, collect(physicalEntity) as peSet\n" +
            "ORDER BY isoform, ptms";

    static final String GET_PROTEIN_NAMES = "MATCH (re:ReferenceEntity{databaseName:'UniProt'})<-[:referenceEntity]-(ewas:PhysicalEntity{speciesName:'Homo sapiens'})\n" +
            "WHERE re.description IS NOT NULL\n" +
            "RETURN DISTINCT \n" +
            "\tCASE WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier ELSE re.identifier END as identifier, \n" +
            "    re.displayName as displayName, \n" +
            "    re.description as description";

    static final String GET_REACTION_PARTICIPANTS = "MATCH p = (rle:ReactionLikeEvent{speciesName:'Homo sapiens'})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nWITH *, relationships(p) as role\nWITH DISTINCT rle.stId as reaction, re.identifier as protein, head(extract(x IN role | type(x))) as role ORDER BY role\nRETURN DISTINCT reaction, protein";

    static final String GET_COMPLEX_PARTICIPANTS = "MATCH p = (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nRETURN DISTINCT c.stId as complex, re.identifier as protein";

}
