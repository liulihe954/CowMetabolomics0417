#load data
data(gse16873.d)
data(demo.paths)
#KEGG view: gene data only
getwd()
i <- 1
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
                   kegg.native = TRUE)
str(pv.out)
head(pv.out$plot.data.gene)
#result PNG file in current directory
#Graphviz view: gene data only
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
                   kegg.native = FALSE, sign.pos = demo.paths$spos[i])
#result PDF file in current directory

#KEGG view: both gene and compound data
sim.cpd.data = sim.mol.data(mol.type="cpd", nmol=3000)
i <- 3
print(demo.paths$sel.paths[i])
pv.out <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix =
                     "gse16873.cpd", keys.align = "y", kegg.native = TRUE, key.pos = demo.paths$kpos1[i])
str(pv.out)
head(pv.out$plot.data.cpd)
length(sim.cpd.data)

#multiple states in one graph
set.seed(10)
sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000,
                              replace = TRUE), ncol = 6)
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y",
                   kegg.native = TRUE, match.data = FALSE, multi.state = TRUE, same.layer = TRUE)
str(pv.out)
head(pv.out$plot.data.cpd)
#result PNG file in current directory
##more examples of pathview usages are shown in the vignette.

?sim.mol.data()


data(cpd.accs)
data(rn.list)
names(rn.list)
cpd.accs[rn.list[[1]][1:4],]
lapply(rn.list[1:4], function(rn) cpd.accs[rn[1:4],])
data(kegg.met)
head(kegg.met)

#kegg.met_test =data.frame(kegg.met)
#kegg.met_test$name  = ""

cpdidmap(in.ids, in.type, out.type)
cpd2kegg(in.ids, in.type)
cpdkegg2name(in.ids, in.type = c("KEGG", "KEGG COMPOUND accession")[1])
cpdname2kegg(in.ids)

