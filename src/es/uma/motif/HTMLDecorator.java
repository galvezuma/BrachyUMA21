package es.uma.motif;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * This program loads a fasta file of proteins, searches for motives and generates a piece of HTML
 * to visualize clearly where are the motives in each sequence.
 * 
 * To make it straightforwardly usable, a CSS piece of HTML is also incorporated. This CSS should
 * be adapted by the developer to match his/her actual needs.
 * 
 */
public class HTMLDecorator {
	
    protected static class ResultPair {
        int pos;
        int qty;
    }
    
    public static void main(String[] args){
    	process("datafiles/all_varieties_DHN_peptide.fasta",
    	                new PerfectMatchArgument[]{
    	                    new PerfectMatchArgument("EKKGIMDKIKEKLPG", "KSegment", 4),
    	                    new PerfectMatchArgument("VDEYGNP", "YSegment", 3),
    	                    new PerfectMatchArgument("EDDGQGR", "InterSegment", 2),
    	                    new PerfectMatchArgument("KKDKKKKKEKK", "NLSSegment", 2),
    	                    new PerfectMatchArgument("DRGLFDKFIGKK", "FSegment", 4), // For Brachypodium
    	                },
    	                new PatternArgument[]{
    	                    //new PatternArgument(Pattern.compile(".K.G..[DE]KIK[DE]K.PG"), "KSegment"),
    	                    //new PatternArgument(Pattern.compile("D[DE][YHF]GNP."), "YSegment"),
    	                    //new PatternArgument(Pattern.compile("LHR[ST]GS{4,6}[SDE][DE]{3}"), "SSegment"),
    	                    new PatternArgument(Pattern.compile("SSSS+"), "SSegment"),
    	                });
    }
    
    public static void process(String filename, PerfectMatchArgument[] pf, PatternArgument[] p) {
        List<Sequence> seqs = loadFasta(filename);
        Map<String, List<Decorator>> ret = new HashMap<>();
        System.out.println("Searching...");
        for(PerfectMatchArgument pfArg: pf) 
            search(seqs, pfArg.aminoAcidSequence, ret, pfArg.name, pfArg.numMismatches);
        for(PatternArgument pArg: p) 
            search(seqs, pArg.pattern, ret, pArg.name);
        seqs.stream().forEach(seq -> {
            int last = seq.data.length();
            String dataDecorated = "";
            if (ret.containsKey(seq.name)) {
                List<Decorator> l = ret.get(seq.name);
                l.sort(null);
                for(int i=l.size()-1; i>=0; i--){
                    if (l.get(i).end <= last) {
                    dataDecorated = 
                            "<span class=\""+l.get(i).classname+"\">"
                            + seq.data.substring(l.get(i).init, l.get(i).end)
                            + "</span>"       
                            + seq.data.substring(l.get(i).end, last)
                            + dataDecorated;
                    last = l.get(i).init;
                    }
                }
            }
            dataDecorated = seq.data.substring(0, last) + dataDecorated;
            System.out.println("<tr><td>&nbsp;</td><td>&nbsp;</td><td class=\"Courier\">"+seq.name+"<br>");
            System.out.println(dataDecorated+"</td></tr>");
        });
        displaySummary(seqs, ret);
    }

    private static void search (List<Sequence> seqs, Pattern r, Map<String, List<Decorator>> ret, String classname) {
        System.out.println(classname);
        seqs.stream().forEach(seq -> {
            int timesFound = 0;
            Matcher m = r.matcher(seq.data);
            if (m.find()){
                List<Decorator> l = ((ret.containsKey(seq.name))? ret.get(seq.name) : new LinkedList<>());
                ret.putIfAbsent(seq.name, l);
                do {
                    l.add(new Decorator(m.start(), m.end(), classname, null));
                    timesFound++;
                } while(m.find());
            }
            System.out.println(seq.name.substring(1)+"\t"+timesFound);
        });
    }
    
    private static void search (List<Sequence> seqs, String pattern, Map<String, List<Decorator>> ret, String classname, int maxBad) {
        System.out.println(classname);
        seqs.stream().forEach(seq -> {
            int timesFound = 0;
            List<Decorator> l = ((ret.containsKey(seq.name))? ret.get(seq.name) : new LinkedList<>());
            ret.putIfAbsent(seq.name, l);
            for(int i=0; i<seq.data.length()-pattern.length()+1; i++){
                String slice = seq.data.substring(i, i+pattern.length());
                boolean[] bol = new boolean[pattern.length()];
                if (match(slice, pattern, bol, maxBad)) {
                    for (int j=0; j<pattern.length(); j++){
                      l.add(new Decorator(i+j, i+j+1, (bol[j])?classname:"aminoDiff", null));  
                    }
                    i += pattern.length() - 1;
                    timesFound++;
                }
            }
            System.out.println(seq.name.substring(1)+"\t"+timesFound);
        });
    }

    private static boolean match(String slice, String pattern, boolean[] bol, int maxBad) {
        int bad = 0;
        for(int i = 0; i < slice.length(); i++){
            bol[i] = (slice.charAt(i) == pattern.charAt(i));
            if (!bol[i] && (++bad) > maxBad) return false;
        }
        return true;
    }

    private static void displaySummary(List<Sequence> seqs, Map<String, List<Decorator>> ret) {
        // Colapse
        boolean continuar = true; // In Spanish because continue is a keyword in Java
        while (continuar){
            continuar = false;
            for(List<Decorator> l: ret.values())
                for(int i=1; i<l.size(); i++) {
                    if(l.get(i).classname.equals("aminoDiff") && l.get(i-1).end == l.get(i).init) l.get(i).classname = l.get(i-1).classname;
                    if(l.get(i-1).classname.equals("aminoDiff") && l.get(i-1).end == l.get(i).init) l.get(i-1).classname = l.get(i).classname;
                    if(l.get(i).classname.equals(l.get(i-1).classname) && l.get(i-1).end == l.get(i).init) {
                        l.get(i-1).end = l.get(i).end;
                        l.remove(i);
                        continuar = true;
                        break;
                    }
                }
        }
        //
        seqs.stream().forEach(seq -> {
            System.out.println(seq.name+"\t"+seq.data.length());
            if (ret.containsKey(seq.name)) {
                List<Decorator> l = ret.get(seq.name);
                for(int i=0; i<l.size(); i++){
                    System.out.println("\t"+l.get(i).classname+"\t"+l.get(i).init+"\t"+l.get(i).end);
                }
            }
        });
    }
    
    public static List<Sequence> loadFasta(String filename) {
        List<Sequence> ret = new ArrayList<>();
        try (BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(filename))))) {
            ret = loadFasta(in);
        } catch (Exception x) {
            x.printStackTrace();
        }
        return ret;
    }
    
    public static List<Sequence> loadFastaGZ(String filename) {
        List<Sequence> ret = new ArrayList<>();
        try (BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename))))) {
            ret = loadFasta(in);
        } catch (Exception x) {
            x.printStackTrace();
        }
        return ret;
    }
    
    public static List<Sequence> loadFasta(BufferedReader in) throws IOException {
        List<Sequence> ret = new ArrayList<>();
        String line;
        StringBuilder sb = new StringBuilder();
        Sequence seq = null;
        while ((line = in.readLine()) != null) {
            if (line.startsWith(">")) {
                if (seq != null) {
                    seq.data = sb.toString();
                    ret.add(seq);
                    sb = new StringBuilder();
                }
                seq = new Sequence(line);
            } else {
                sb.append(line);
                //seq.data += line;
            }
        }
        if (seq != null) {
            seq.data = sb.toString();
            ret.add(seq);
            sb = null;
        }
        return ret;
    }
    
}


