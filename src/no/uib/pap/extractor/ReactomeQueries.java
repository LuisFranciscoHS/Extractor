package no.uib.pap.extractor;

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
    
    public static final String GET_ALL_PATHWAYS = "MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nWITH DISTINCT  pe, re\nOPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)\nWITH DISTINCT \n  pe,\n  re.identifier AS protein,\n  re.variantIdentifier AS isoform,\n  tm.coordinate as coordinate, \n  mod.identifier as type ORDER BY type, coordinate\nWITH DISTINCT   \n  pe,\n  protein,\n  isoform,\n  COLLECT(CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE \"null\" END + \":\" + type) AS ptms \n WITH DISTINCT   \n  pe,\n  protein,\n  COLLECT(CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END + ptms) as proteoform\n  MATCH (p:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName: \"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe)\nRETURN DISTINCT p.stId AS stId, count(DISTINCT rle) as numReactionsTotal, count(DISTINCT protein) as numEntitiesTotal, count(DISTINCT proteoform) as numProteoformsTotal\n  ORDER BY stId";

    public static final String GET_COUNT_ALL_PATHWAYS = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"}) RETURN count(DISTINCT p) as count";

    public static final String GET_MAP_REACTIONS_TO_PATHWAYS = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName: \"Homo sapiens\"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:\"UniProt\"})\nRETURN DISTINCT p.stId AS pathway, rle.stId as reaction";
    
    public static final String GEP_MAP_PATHWAYS_TO_TOPLEVELPATHWAYS = "MATCH (tlp:TopLevelPathway{speciesName:'Homo sapiens'})-[:hasEvent*]->(p:Pathway{speciesName:'Homo sapiens'})\n" + 
    		"RETURN DISTINCT p.stId AS pathway, tlp.stId as topLevelPathway";

	public static final String GET_MAP_PROTEOFORMS_TO_PHYSICALENTITIES = "MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})\nWITH DISTINCT pe, re\nOPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)\nWITH DISTINCT pe.stId AS physicalEntity,\n                re.identifier AS protein,\n                re.variantIdentifier AS isoform,\n                tm.coordinate as coordinate, \n                mod.identifier as type ORDER BY type, coordinate\nWITH DISTINCT physicalEntity,\n                CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END as isoform,\n                COLLECT(type + \":\" + CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE \"null\" END) AS ptms\n                RETURN isoform, ptms, collect(physicalEntity) as peSet\n                ORDER BY isoform, ptms";

}
