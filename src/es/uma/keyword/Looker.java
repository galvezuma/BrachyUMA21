package es.uma.keyword;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;
import es.uma.motif.Sequence;
import es.uma.motif.ResultPair;

/**
 *
 * @author galvez
 */
public class Looker {
    private static boolean firstLine = true;
    public static StringBuilder str= new StringBuilder(750_000_000);
    public static HashMap<String, ResultPair> hm = new HashMap<>();

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

    public static void loadChromosome(String name) {
        str.delete(0, Integer.MAX_VALUE);
        try (BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File("G:\\Data_JIC\\IWGSC\\WGAv1.0\\Chromosomes\\chr"+name+".fa"))))) {
            in.lines().forEach(l -> {
                if (firstLine) {
                    firstLine = false;
                } else {
                    str.append(l);
                }
            });
        } catch (Exception x) {
            System.out.println(Runtime.getRuntime().maxMemory());
            x.printStackTrace();
        }
    }
 
}

