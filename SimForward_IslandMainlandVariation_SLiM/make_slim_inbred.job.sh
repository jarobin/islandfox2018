# Script to make SLiM job script for Inbreeding model

# Set number of genes
g=${1}

# Set number of burn-in generations to build up variation in the ancestral population
t=${2}

# Set per-site mutation rate
u=${3}

# Set long-term population size on the island (e.g. 50, 100, 200, 500, 1000)
n=${4}

# Set dominance coefficient for deleterious mutations 
# (e.g. h=0.0 for fully recessive, h=0.5 for additive)
h=${5}

# Specify model
m=inbred


cat > slim_${m}_${g}genes_${t}burn_n${n}_h${h}.job << EOM
// Inbreeding Model

// INITIALIZE
initialize() {
// For inbreeding, must track pedigrees
initializeSLiMOptions(keepPedigrees=T);

initializeMutationRate(${u}); 

// Distribution of selection coefficients for nonsynonymous mutations 
// from Kim et al. 2017, DOI:10.1534/genetics.116.197145
initializeMutationType("m1", ${h}, "g", -0.01314833, 0.186); 

// Neutral mutations
initializeMutationType("m2", 0.5, "f", 0.0);

// Synonymous:nonsynonymous mutation ratio of 1:2.31 from Huber et al. 2017,
// DOI:10.1073/pnas.1619508114
initializeGenomicElementType("g1", c(m1,m2), c(2.31,1.0));

// Make "genes", each 1 kb in length
for (i in 1:${g}){
	initializeGenomicElement(g1, ((i-1)*1000)+(i-1), (i*1000)+(i-2) );
}

// Number of genes per chromosome for 2000 genes, assuming genes are evenly distributed
// and chromosome lengths correspond to relative lengths in the dog genome
// Note: these numbers must be changed if the number of genes is not 2000
gene_nums=c(111,78,83,80,81,70,73,68,55,63,68,66,57,55,58,54,58,51,49,53,46,56,48,43,47,35,42,37,38,37,36,35,29,38,24,28,28,22);

// Set recombination rates:
// Recombination rate in between "chrosomomes" is 50% (i.e. unlinked)
// Recombination rate within genes is 0
// Recombination rate in between genes is 1e-3, the effective recombination rate for 
// 100 kb noncoding sequence with a per-site recombination rate of 1e-8 in between 
// each gene
rates=NULL;
for (i in 1:(size(gene_nums)-1)){
	rates=c(rates, 0, rep(c(1e-3, 0), gene_nums[i-1]-1), 0.5);
}
rates=c(rates, 0, rep(c(1e-3, 0), gene_nums[size(gene_nums)-1]-1));

ends=NULL;
for (i in 1:${g}){
	ends=c(ends, (i*1000)+(i-2), (i*1000)+(i-1));
}
ends=ends[0:(size(ends)-2)];

initializeRecombinationRate(rates, ends);

}


// SIMULATE
// Begin the simulation with burn-in period
1 { 
	sim.addSubpop("p1", 10000); 
} 

// After burn-in, create island population
${t} {
    sim.addSubpopSplit("p2", ${n}, p1);
}

// Inbreeding: closely related individuals (relatedness >= 0.25), are twice as likely to 
// mate with one another, but no selfing
$((${t} + 10000 - 100)): mateChoice(p2) {
	if (runif(1) < 0.67) {
			weights * ifelse(individual.relatedness(sourceSubpop.individuals) >= 0.25 & individual.relatedness(sourceSubpop.individuals) != 1.0 , 1.0 , 0.00001);
	}
	else {
		weights;
	}
}

// Keep track of generation number in log file every 1,000 generations (overwrites itself)
1:$((${t} + 10000)) late() {
	if (sim.generation % 1000 == 0){
		writeFile("slim_${m}_${g}genes_${t}burn_n${n}_h${h}.job.gen", paste(sim.generation));
	}
}

// After 10,000 generations the simulation ends and output is written to stdout
$((${t} + 10000)) late() {
  // Output population sizes for easy reference
  cat("#p1: " + size(p1.individuals) + " individuals; p2: " + size(p2.individuals) + " individuals.\n");
  
  // File header:
  // Mutation id,
  // Mutation type,
  // Selection coefficient,
  // Age of mutation in generations,
  // Subpopulation it arose in,
  // Number of heterozygote derived in p1,
  // Number of homozygote derived in p1,
  // Number of heterozygote derived in p2,
  // Number of homozygote derived in p2
  // Note: these are genotype counts, not allele counts
  cat("mutid,type,s,age,subpop,p1numhet,p1numhom,p2numhet,p2numhom\n");
  
  // Collect stats for each mutation in the simulation
  for (mut in sim.mutations){
    id = mut.id;
    s = mut.selectionCoeff;
    subpop = mut.subpopID;
    age = sim.generation - mut.originGeneration;
    type = mut.mutationType;
    if (type == m1){
      type = "m1";
    } else if (type == m2){
      type = "m2";
    }

    // Initialize genotype counts
    p1numhet = 0;
    p1numhom = 0;
    p2numhet = 0;
    p2numhom = 0;
    
    // Count number of derived heterozygotes and homozygotes in p1
    for (p1i in p1.individuals){
      gt = sum(c(p1i.genomes[1].containsMutations(mut), p1i.genomes[0].containsMutations(mut)));
      if (gt == 1){
        p1numhet = p1numhet + 1;
      } else if (gt == 2){
        p1numhom = p1numhom + 1;
      }
    }
    
    // Count number of derived heterozygotes and homozygotes in p2
    for (p2i in p2.individuals){
      gt = sum(c(p2i.genomes[1].containsMutations(mut), p2i.genomes[0].containsMutations(mut)));
      if (gt == 1){
        p2numhet = p2numhet + 1;
      } else if (gt == 2){
        p2numhom = p2numhom + 1;
      }
    }

    // Print results
    cat(id + ",");
    cat(type);
    cat("," + s + "," + age + "," + subpop + "," + p1numhet + "," + p1numhom + "," + p2numhet + "," + p2numhom + "\n");
  }
}

EOM
