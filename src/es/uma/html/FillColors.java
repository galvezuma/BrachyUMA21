package es.uma.html;

import nw.NWAlign;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class FillColors {
    protected static final String[] COLORS = { 
        "#66FFFF", // Bdhn10
        "#FFCCFF", // Bdhn3
        "#9AFED6", // Bdhn9
        "#99FF99", // Bdhn8
        "#CCFF99", // Bdhn7
        "#CCCCFF", // Bdhn2
        "#FFFF99", // Bdhn6
        "#FFCCCC", // Bdhn4
        "#FFCC99", // Bdhn5
        "#66CCFF", // Bdhn1
    };
    // This array indicates at position i-th where is stored the colour of the Bdhn i-th in the previous array.
    // Indexes are consistent with numbers of dehydrins and, therefore, they begin with 1 and not with 0.
    protected static final int[] POSITIONS = {10, 6, 2, 8, 9, 7, 5, 4, 3, 1};
    private static final String basicColor = "#8B0000"; // DarkRed
    private static final String sequenceMissingColor = "#FF0000"; // Red
    
    public static void iterateChromosomes(int chrNum, List<Variety> variedades) {
        int indexColor = -1;
        for(Variety vModelo: variedades) {
            for(Gene g: vModelo.genes){
                if (g.getChrNumber() != chrNum) continue;
                if (g.color != null && ! g.color.equals(basicColor)) continue;
                if (g.protein.equals("")) { g.color = sequenceMissingColor; continue; }
                g.color = COLORS[++indexColor];
                for(Variety v: variedades)
                    if (v != vModelo) iterateVariety(v, chrNum, g, indexColor);
            }
        }
        System.out.println(indexColor);
    }

    private static void iterateVariety(Variety v, int chr, Gene g, int indexColor) {
        // Gene bestMatch = null;
        double optimalScore = (double)needlemanWunsch(g.protein, g).score;
        // double bestPercentage = 0.0;
        
        ExecutorService executorService = Executors.newFixedThreadPool(8);
        List<Future<NWResult>> execs = new ArrayList<>();
        for(Gene g2: v.genes) if (g2.getChrNumber() == chr && g2.color == null && !g.protein.equals("")) 
            execs.add(executorService.submit(() -> needlemanWunsch(g.protein, g2)));
        
        for(Future<NWResult> f: execs) try {
            NWResult result = f.get();
            if (result.score/optimalScore > 0.95) {
                result.gene.color = COLORS[indexColor];
                result.gene.info += "\t" + result.score/optimalScore;
//                if (bestMatch == null || Math.abs(result.gene.init-g.init) < Math.abs(bestMatch.init-g.init)) {
//                    bestMatch = result.gene;
//                    bestPercentage = result.score/optimalScore;
//                } 
            }
        } catch (Exception x) { x.printStackTrace(); }
        executorService.shutdownNow();
//        if (bestMatch != null) {
//            bestMatch.color = COLORS[indexColor];
//            bestMatch.info += "\t" + bestPercentage;
//        }
    }
    
    /* 
     * Once coloured the similar genes of a particular chromosome in all the varieties, we address the problem
     * of finding out if there are similar genes in unknown chromosomes and drawing them in the most likely
     * chromosome, even if this should be done at position 0.
     * To do this, for each colour of the just coloured chromosome (1, 2, etc.), let's iterate over the first 
     * gene associated to each colour and let's look for it in the set of unknown chromosomes of those varieties
     * without such a colour in the chromosome just coloured. We will change the chromosome of these genes in order 
     * to associate them to the just coloured chromosome.
    */
    public static void iterateUnknownChromosomes(int chrNum, List<Variety> variedades) {
        Set<String> usedColors = new HashSet<>();
        for(Variety vModelo: variedades) {
            for(Gene g: vModelo.genes){
                if (g.getChrNumber() != chrNum) continue;
                if (g.color == null) continue;
                if (g.color.equals(sequenceMissingColor)) continue;
                if (usedColors.contains(g.color)) continue;
                usedColors.add(g.color);
                double optimalScore = (double)needlemanWunsch(g.protein, g).score;
                System.out.println("Processing "+g.color+" with score "+optimalScore);
                // Reassign genes from Unknown chromosomes
                ExecutorService executorService = Executors.newFixedThreadPool(8);
                List<Future<NWResult>> execs = new ArrayList<>();
                for(Variety variety: variedades) {
                    if (hasColorInChromosome(variety, g.color, chrNum)) continue;
                    for(Gene g2: variety.genes) if (g2.getChrNumber() == -1) {
                        execs.add(executorService.submit(() -> needlemanWunsch(g.protein, g2)));
                        System.out.println("Trying "+g2.name);
                    }
                }
                for(Future<NWResult> f: execs) try {
                    NWResult result = f.get();
                    if (result.score/optimalScore > 0.95) {
                        result.gene.color = g.color;
                        result.gene.info += "\t" + result.score/optimalScore;
                        result.gene.chr = ""+chrNum;
                        result.gene.init = 0;
                        result.gene.shape = "CIRCLE";
                        System.out.println("Altered "+result.gene.name);
                    }
                } catch (Exception x) { x.printStackTrace(); }
                executorService.shutdownNow();
                //
            }
        }
    }

    private static boolean hasColorInChromosome(Variety variety, String color, int chrNum) {
        for(Gene g: variety.genes) if (g.getChrNumber() == chrNum) {
            if (color.equals(g.color)) return true;
        }
        return false;
    }

    protected static NWResult needlemanWunsch(String q, Gene s) {
        int gapOpen=-11, gapExtn=-1;
        NWResult ret = new NWResult();
        ret.gene = s;
        ret.score = NWAlign.NeedlemanWunsch(q, s.protein, gapOpen, gapExtn);
        return ret;
    }

    static void iterateVarietiesAgainstFirst(List<Variety> variedades) {
        int indexColor = -1;
        // A different colour is assigned to each gene of the model (first variety, usually B. distachyon).
        Variety vModelo = variedades.get(0);
        for(Gene g: vModelo.genes){
            g.color = COLORS[++indexColor];
        }
        // For any other gene of any other variety, we assign to it the colour
        // of the most similar gene in the model.
        // In addition, if the chromosome of such a variety gene is unknown, we assign
        // the same chromosome of the gene selected in the model.
        for(int i=1; i<variedades.size(); i++) {
            for(Gene g: variedades.get(i).genes){
                Gene geneModelo = mostSimilarGene(g, vModelo);
                g.color = geneModelo.color;
                if (g.getChrNumber() == -1) {
                    g.chr = geneModelo.chr;
                    double scale = (double)(variedades.get(i).chrLength[geneModelo.getChrNumber()-1])/(double)(vModelo.chrLength[geneModelo.getChrNumber()-1]);
                    g.init = (int)(geneModelo.init*scale); // Estimate position;
                    g.shape = "CIRCLE";
                    System.out.println("Altered "+g.name+" "+g.init);
                }
            }
        }
    }

    // This function returns the most similar gene in the model.
    // The score given by N/W is used to find out such a most similar gene.
    // This is done in parallel.
    protected static Gene mostSimilarGene(final Gene g, final Variety vModelo) {
        // Prepare parallel execution
        ExecutorService executorService = Executors.newFixedThreadPool(8);
        List<Future<NWResult>> execs = new ArrayList<>();
                
        execs.add(executorService.submit(() -> needlemanWunsch(g.protein, vModelo.genes[0])));
        assert(! execs.isEmpty());
        // Next is the same than for(int i=1; i<vModelo.genes.length; i++)
        // but creates immutable objects i that can be used inside
        // a lambda expression.
        for(Integer i: IntStream.rangeClosed(1, vModelo.genes.length-1).boxed().collect(Collectors.toList())) {
            execs.add(executorService.submit(() -> needlemanWunsch(g.protein, vModelo.genes[i])));
        }
        NWResult bestResult = null;
        for(Future<NWResult> f: execs) try {
            if (bestResult == null) { bestResult = f.get(); }
            else if (f.get().score > bestResult.score) {
                bestResult = f.get();
            }
        } catch (Exception x) { x.printStackTrace(); }
        executorService.shutdownNow();
        return bestResult.gene;
    }
    
    protected static class NWResult {
        int score;
        Gene gene;
    }
}
