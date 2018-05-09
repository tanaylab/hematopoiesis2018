require("devtools");
load_all("metacell/tgconfig")
load_all("metacell/tgstat/")
load_all("metacell/tgutil/")
load_all("metacell/scrdb/")
set_param("use_tgs_cor", F, package = "scrdb")

source("scripts/sc_functions.r")
set.seed(2754)

index_fn = "annotations/batch_stats.txt"
index = read.delim(index_fn, stringsAsFactors = F)
batch_stats = read.delim(index_fn, stringsAsFactors = F, row.names = 1)

genes = rownames(read.delim("output/umi.tab/AB1154.txt", row.names = 1))
modules = read.delim("annotations/modules.txt", row.names = 1)
cc_genes = rownames(modules)[ modules$annotation == "CC"]
ribo_genes = union(rownames(modules)[ modules$annotation == "Ribo"], grep("Rpl|Rps", genes, value = T))
bad_genes = c("Atpase6", cc_genes, ribo_genes)

dir.create("saved_work")
dir.create("figures")
dir.create("supp_figures")
dir.create("temp")

wells_cells = "annotations/metadata.txt"
wells = read.delim(wells_cells, stringsAsFactors = F, row.names = 1, skip=40)
#plate_facs = read.delim("annotations/facs_coordinates.txt", stringsAsFactors = F, row.names = 2)[,-1]
load("gates/alternative_gates.Rda")
lin_gate = gate_facs(wells, "Lineage", "cKit", log = "xy", rect = T, polygon=lin_gate$polygon)
X = wells[ rownames(wells), "Subject_ID"] == "FP76"
lin_gate$gate = (X & wells$Lineage < 200) | (!X & wells$Lineage < 1000)
names(lin_gate$gate) = rownames(wells)
ckit_gate = gate_facs(wells, "Lineage", "cKit", log = "xy", rect = T, polygon=ckit_gate$polygon)
names(ckit_gate$gate) = rownames(wells)

load("gates/tier3_gates.Rda")
lsk = gate_facs(wells, "Sca1", "cKit", log = "xy", rect = T, polygon=lsk$polygon)
lk = gate_facs(wells, "Sca1", "cKit", log = "xy", rect = T, polygon=lk$polygon)
lt = gate_facs(wells, "CD34", "Flt3", log = "xy", rect = T, polygon=lt$polygon)
st = gate_facs(wells, "CD34", "Flt3", log = "xy", rect = T, polygon=st$polygon)
mpp = gate_facs(wells, "CD34", "Flt3", log = "xy", rect = T, polygon=mpp$polygon)
mep = gate_facs(wells, "CD34", "FcgR", log = "xy", rect = T, polygon=mep$polygon)
cmp = gate_facs(wells, "CD34", "FcgR", log = "xy", rect = T, polygon=cmp$polygon)
gmp = gate_facs(wells, "CD34", "FcgR", log = "xy", rect = T, polygon=gmp$polygon)
clp = gate_facs(wells, "IL7Ra", "cKit", log = "xy", rect = T, polygon=clp$polygon)

tier3_gates = data.frame(well = rownames(wells), lsk = lsk$gate, lt = lt$gate, st = st$gate,
                   mpp = mpp$gate, lk = lk$gate, mep = mep$gate,
                   cmp = cmp$gate, gmp = gmp$gate, clp = clp$gate)
rownames(tier3_gates) = tier3_gates$well

tier3_gates$gate = "other"
tier3_gates$gate[ tier3_gates$lsk & !is.na(tier3_gates$lsk)] = 
		  colnames(tier3_gates)[3:5][apply(cbind(tier3_gates[ tier3_gates$lsk & !is.na(tier3_gates$lsk), 3:5], T), 1, which.max)]
tier3_gates$gate[ tier3_gates$lk & !is.na(tier3_gates$lk)]   = 
		  colnames(tier3_gates)[7:9][apply(cbind(tier3_gates[ tier3_gates$lk & !is.na(tier3_gates$lk),  7:9],T) , 1, which.max)]
tier3_gates$gate[ tier3_gates$clp] = "clp"

load("gates/tier6_gates.Rda")
lt = gate_facs(wells, "CD34", "Flt3", log = "xy", rect = T, polygon=lt$polygon)
st = gate_facs(wells, "CD34", "Flt3", log = "xy", rect = T, polygon=st$polygon)
mpp = gate_facs(wells, "CD34", "Flt3", log = "xy", rect = T, polygon=mpp$polygon)
hsc = gate_facs(wells, "CD48", "CD150", log = "xy", rect = T, polygon=hsc$polygon)
mpp1 = gate_facs(wells, "CD48", "CD150", log = "xy", rect = T, polygon=mpp1$polygon)
mpp2 = gate_facs(wells, "CD48", "CD150", log = "xy", rect = T, polygon=mpp2$polygon)
mpp3 = gate_facs(wells, "CD48", "CD150", log = "xy", rect = T, polygon=mpp3$polygon)
mpp4 = gate_facs(wells, "CD48", "CD150", log = "xy", rect = T, polygon=mpp4$polygon)

tier6_gates = data.frame(well = rownames(wells), lt = lt$gate, st = st$gate,
                   mpp = mpp$gate, hsc = hsc$gate, mpp1 = mpp1$gate, mpp2 = mpp2$gate,
                   mpp3 = mpp3$gate, mpp4 = mpp4$gate)
rownames(tier6_gates) = tier6_gates$well
tier6_gates$norm_gate = "other"
tier6_gates$norm_gate = colnames(tier6_gates)[2:4] [ apply(cbind(tier6_gates[, 2:4], T), 1, which.max)]

tier6_gates$new_gate = "other"
tier6_gates$new_gate[ tier6_gates$lt & !is.na(tier6_gates$lt)] = ifelse(tier6_gates[ tier6_gates$lt & !is.na(tier6_gates$lt), "hsc"], "hsc", "other")
tier6_gates$new_gate[ tier6_gates$mpp & !is.na(tier6_gates$mpp)] = ifelse(tier6_gates[ tier6_gates$mpp & !is.na(tier6_gates$mpp), "mpp4"], "mpp4", "other")
tier6_gates$new_gate[ tier6_gates$st & !is.na(tier6_gates$st)] = 
		      colnames(tier6_gates)[6:8] [ apply(cbind(tier6_gates[tier6_gates$st & !is.na(tier6_gates$st), 6:8], T), 1, which.max)]


#gates[is.na(gates)] = "other"

load("gates/mye_gates.Rda")
flt3_gate    = gate_facs(wells, "Csf1r",    "Flt3", log = "xy", rect = T, polygon=flt3_gate$polygon)
dp_gate      = gate_facs(wells, "Csf1r",    "Flt3", log = "xy", rect = T, polygon=dp_gate$polygon)
csf1r_gate   = gate_facs(wells, "Csf1r",    "Flt3", log = "xy", rect = T, polygon=csf1r_gate$polygon)
mdp_gate     = gate_facs(wells, "CD11c",    "cKit", log = "xy", rect = T, polygon=mdp_gate$polygon)
cdp_gate     = gate_facs(wells, "CD11c",    "cKit", log = "xy", rect = T, polygon=cdp_gate$polygon)
predc_gate   = gate_facs(wells, "CD11c",    "cKit", log = "xy", rect = T, polygon=predc_gate$polygon)
cmop_gate    = gate_facs(wells, "Ly6C",     "cKit", log = "xy", rect = T, polygon=cmop_gate$polygon)
pdc_gate     = gate_facs(wells, "SiglecH",  "Flt3", log = "xy", rect = T, polygon=pdc_gate$polygon)

mye_gates = data.frame(well = rownames(wells), flt3 = flt3_gate$gate, dp = dp_gate$gate,
            csf1r = csf1r_gate$gate, mdp = mdp_gate$gate, cdp = cdp_gate$gate, predc = predc_gate$gate,
            cmop = cmop_gate$gate, pdc = pdc_gate$gate, clp = clp$gate)
rownames(mye_gates) = mye_gates$well
mye_gates$gate = "other"

ind = which(!is.na(mye_gates$dp))
mye_gates[ ind, "gate"] = with(mye_gates[ind,],
	ifelse(flt3 & pdc,     "pDC",
	ifelse(csf1r & cmop,   "cMoP", 
	ifelse(dp & mdp,       "MDP",
	ifelse(dp & cdp,       "CDP",
	ifelse(flt3 & cdp,     "Csf1r- CDP",
	ifelse(dp & predc,     "preDC", "other")))))))



