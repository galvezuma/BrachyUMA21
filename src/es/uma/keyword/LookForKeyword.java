package es.uma.keyword;



import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import es.uma.motif.Sequence;
import es.uma.html.Gene;

/**
 * This program extract dehydrins from annotation files taken from Phytozome.
 * 
 * It displays information for each Brachypodium distachyon variety:
 * - Name of the variety.
 * - Total number of dehydrins found.
 * - Length of each chromosome.
 * 
 * For each dehydrin, it displays:
 * - Data taken from the annotation file.
 * - Gene sequence.
 * - Data taken from the gff3 file.
 * - Peptide sequence.
 * The result should copied and pasted into "dehydrins.sal" for further processing.
 * 
 * @author galvez
 */
public class LookForKeyword {
    private static final String DIR_NAME = "genomes/Phytozome/PhytozomeV13";
    private static final String INPUT_DATA_FILE = "genomes/Order.txt";
    public static void main(String[] args) throws Exception {
        Path[] filesInOrder = getFilesInOrder();
        for(Path file: filesInOrder)
            process(file); 
    }
    
    private static void process(Path file) throws Exception {
        HashMap<String, String> cjt = new HashMap<>();
        String geneName = null;
        int counter=0;
        System.out.println(file.getFileName().toString());
        try {
            for (String line : Files.readAllLines(file, StandardCharsets.UTF_8))
            	// Change these to look for other families of genes
                if (line.toUpperCase().contains("PF00257") || 
                    line.toUpperCase().contains("DEHYDRIN")) {
                    String transcriptName = line.split("\t")[2];
                    geneName = removeTranscriptNumber(transcriptName);
                    if(cjt.containsKey(geneName)) {
                        // System.err.println("Repeated: "+transcriptName);
                    } else {
                        counter++;
                        String shape = line.toUpperCase().contains("PF00257")?"RECT":"TRIANGLE";
                        // System.out.println("OK:\t"+transcriptName);
                        cjt.put(geneName, line+"\t"+shape);
                    }
                }
        } catch (IOException ex) { ex.printStackTrace(); }
        System.out.println("Total:\t"+counter);
        Map<String, Gene> genesAll = loadGff3(file.getParent(), (geneName.endsWith("m")?"mRNA":"gene")); 
        displayChrLengths(genesAll);
        List<Sequence> seqsAll = loadFastaPrimaryTranscript(file.getParent());
        List<Sequence> genesGenome = loadFastaGenesGenome(file.getParent());
        displaySorted(cjt, genesAll, seqsAll, genesGenome);
    }

    private static String removeTranscriptNumber(String transcriptName) {
        if (! transcriptName.contains(".")) return transcriptName;
        String[] pieces = transcriptName.split("\\.");
        String geneName = pieces[0];
        for(int i=1; i<pieces.length-1; i++){
            geneName += "." + pieces[i];
        }
        return geneName;
    }
    
    public static List<Sequence> loadFastaPrimaryTranscript(Path parent) {
        return loadFasta(parent, ".protein_primaryTranscriptOnly.fa.gz");
    }
    public static List<Sequence> loadFastaGenesGenome(Path parent) {
        return loadFasta(parent, "genes.fasta");
    }
    private static List<Sequence> loadFasta(Path parent, String end) {
        List<Sequence> ret = null;
        try {
            File[] fileList = parent.toFile().listFiles((dir, name) -> name.endsWith(end));
            BufferedReader in = null;
            if (fileList[0].getAbsolutePath().endsWith(".gz"))
                in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileList[0]))));
            else 
                in = new BufferedReader(new FileReader(fileList[0]));
            ret = Looker.loadFasta(in);
            in.close();
        } catch (IOException ex) { ex.printStackTrace(); }
        return ret;
    }
    
    private static Map<String, Gene> loadGff3(Path parent, String toSearch) {
        HashMap<String, Gene> cjto = new HashMap<>();
        try {
            File[] fileList = parent.toFile().listFiles((dir, name) -> name.endsWith(".gene.gff3.gz"));
            BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fileList[0]))));
            in.lines().forEach(line -> {
                String[] pieces = ((String)line).split("\t");
                if (pieces.length >= 9 && pieces[2].equals(toSearch)){
                    int init = Integer.parseInt(pieces[3]);
                    String name = pieces[8].split(";")[1];
                    name = name.substring(5);
                    String chr = pieces[0]; //.substring(pieces[0].indexOf("_")+1);
                    cjto.put(name, new Gene(name, init, chr, line));
                }
            });
            in.close();
        } catch (IOException ex) { ex.printStackTrace(); }
        return cjto;
    }

    private static void displaySorted(Map<String, String> set, Map<String, Gene> genesAll, List<Sequence> seqsAll, List<Sequence> genesGenome) throws Exception {
        ArrayList<Gene> list = new ArrayList<>();
        set.keySet().stream().forEach(name -> { 
            Gene gen = genesAll.get(name);
            if (gen == null) System.err.println("MISSING IN GFF3:\t"+name); 
            else {
                gen.info = set.get(name);
                list.add(gen);
            } 
        } );
        for(Gene g: list) {
            g.protein = locate(g.name, seqsAll).data;
            g.gene = locate(g.name, genesGenome).data;
        }
        Collections.sort(list);
        list.stream().forEach(System.out::println);
    }

    private static void displayChrLengths(Map<String, Gene> genesAll) {
        int[] lengths = new int[5];
        for(Gene g: genesAll.values()){
            int chr = g.getChrNumber();
            if (chr != -1 && lengths[chr-1]<g.init) lengths[chr-1] = g.init;
        }
        for(int i=0; i<5; i++)
            System.out.println((i+1)+"\t"+lengths[i]);
    }

    private static Sequence locate(String name, List<Sequence> seqsAll) throws Exception {
        Sequence ret = null;
        for(Sequence seq: seqsAll){
            if (seq.name.startsWith(">"+name)) { // contains(name)) {
                if (ret == null) ret = seq;
                else throw new Exception("DUPLICATED PROTEIN:\t"+name);
            }
        }
        if (ret == null) ret = Sequence.getEmptySequence(); // Some sequences are missing in the fasta file
        return ret;
    }

    /*
    Return the annotation files (with extension .txt) in the order given by the INPUT_DATA_FILE.
    */
    public static Path[] getFilesInOrder() {
        try (BufferedReader in = Files.newBufferedReader(Paths.get(INPUT_DATA_FILE));){
            List<String> order = new ArrayList<>();
            String line;
            while ((line = in.readLine()) != null)
                if (! line.startsWith("#"))
                    order.add(line.split("\t")[0]);
            Path[] ret = new Path[order.size()];
            Files.walkFileTree(Paths.get(DIR_NAME), new SimpleFileVisitor<Path>(){
                @Override
                public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                    String fileName = file.toAbsolutePath().toString();
                    if (fileName.endsWith(".txt")) {
                        int pos = search(order, fileName);
                        if (pos >= 0)
                            ret[pos] = file;
                    }
                    return FileVisitResult.CONTINUE;
                }
            });
            for(int i=0; i<ret.length; i++)
                if (ret[i] == null)
                    throw new Exception("Problem looking for: "+order.get(i));
            System.out.println("Working with "+ret.length+" files.");
            return ret;
        } catch (Exception x) { x.printStackTrace(); }
        return null;
    }

    // Name of ecotype should be part of the filename's path
    //
    private static int search(List<String> order, String fileName) {
        String[] subdirs = fileName.split(Pattern.quote(FileSystems.getDefault().getSeparator()));
        for(String subdir: subdirs){
            if (order.contains(subdir))
                return order.indexOf(subdir);
        }
        return -1;
    }
    
}



