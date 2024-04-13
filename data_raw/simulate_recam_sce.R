## Code to simulate and normalize example data

library(RESCUE)

## Set seed for reproducibility
set.seed=24

## Load present
data("RecAM_params")

## Adjust parameters
RecAM_params <- updateRescueParams(
    paramObj = RecAM_params,
    paramValues = list(
        twoGroupDesign = T,
        maxCellsPerSamp=150,
        minCellsPerSamp=50,
        propDE=.2,
        deLogFC=1
    )
)

## Simulate data
RecAM_sim_sce <- simRescueData(RecAM_params)

## Remove genes that have more than 75% 0's
idx_rm=which(rowSums(assay(RecAM_sim_sce, "counts")==0)>(.75*ncol(RecAM_sim_sce)))
RecAM_sim_sce=RecAM_sim_sce[-idx_rm,]

assay(RecAM_sim_sce, "normcounts")<-sctransform::vst(counts(RecAM_sim_sce))$y

usethis::use_data(RecAM_sim_sce, overwrite = TRUE)



RecAM_sim_pb=rebelAggregate(object=RecAM_sim_sce, countsAssay="counts", sampleVariable="sampleID")

assay(RecAM_sim_pb,"normcounts")=DESeq2::varianceStabilizingTransformation(assay(RecAM_sim_pb,"counts"),
                                                                           blind = F)

usethis::use_data(RecAM_sim_pb, overwrite = TRUE)


