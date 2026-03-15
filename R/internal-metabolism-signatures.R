# Curated metabolism related signatures used by calculate_metabolism_score().
# This object is internal to the package and intentionally not exported.
signature_metabolism <- list(
  Glycolysis = c("HK1", "HK2", "PFKM", "ALDOA", "ENO1", "PKM", "LDHA"),
  Oxidative_Phosphorylation = c("NDUFA1", "NDUFB8", "UQCRC1", "COX4I1", "ATP5F1A", "ATP5F1B"),
  TCA_Cycle = c("CS", "ACO2", "IDH3A", "OGDH", "SUCLG1", "SDHA", "FH", "MDH2"),
  Lipid_Metabolism = c("FASN", "ACACA", "SCD", "CPT1A", "ACOX1", "HMGCR")
)
