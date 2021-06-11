package es.uma.html;

import java.util.Objects;

public class Gene implements Comparable<Object> {
    public String name;
    public String chr;
    public int init;
    public String info="";
    public String color; // DarkRed
    public String protein;
    public String shape;
    public String gff3EntryGene;
    public String gene;
    public String dhnCode; // Code of the reference Bdhn
    
    public Gene(String name, int init, String chr, String gff3){
        this.name = name;
        this.init = init;
        this.chr = chr;
        this.gff3EntryGene = gff3;
        // System.out.println(name+"\t"+chr+"\t"+init);
    }

    public Gene(String line) {
        try {
            String[] pieces = line.split("\t");
            this.name = pieces[0];
            this.chr = pieces[1];
            this.init = Integer.parseInt(pieces[2]);
            this.shape = pieces[pieces.length - 1];
            this.info = line.replace("\t"+this.shape, "");
        } catch (Exception x) { System.err.println(line); x.printStackTrace(); System.exit(0); }
    }
    
    @Override
    public boolean equals(Object o) {
        // self check
        if (this == o) return true;
        // null check
        if (o == null) return false;
        // type check and cast
        if (getClass() != o.getClass()) return false;
        Gene gene = (Gene) o;
        // field comparison
        return Objects.equals(name, gene.name);
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 53 * hash + Objects.hashCode(this.name);
        return hash;
    }

    @Override
    public int compareTo(Object o) {
        if ((o == null) || (! (o instanceof Gene))) return -1;
        else {
            Gene g = (Gene) o;
            if (this.equals(o)) return 0;
            else if (this.chr.compareTo(g.chr) != 0) return this.chr.compareTo(g.chr);
            else if (this.init < g.init) return -1;
            else return +1;
        }
    }
    
    @Override
    public String toString(){
        return name+"\t"+chr+"\t"+init+"\t"+info+"\n"+protein+"\n"+gff3EntryGene+"\n"+gene;
    }

    public int getChrNumber() {
        int ret=0;
        for(int i=0; i<chr.length(); i++){
            if (Character.isDigit(chr.charAt(i))){ 
                ret = ret*10+Character.getNumericValue(chr.charAt(i));
                if (ret>100) break;
            }
        }
        if (ret<1 || ret>5) return -1;
        return ret;
    }
    
    public String getColor() {
        return (color==null)?"#8B0000":color;
    }

    public String getHTML() {
        StringBuilder stb = new StringBuilder();
        stb.append("<table border='1'><tr><td>Chr.</td><td>").append(chr).append("</td></tr>");
        String[] pieces = this.info.split("\t");
        int[] data = new int[]{0, 2, 7, 8, 9, 12, 15, 18, 21, 22};
        String[] info = new String[]{
            "Gene",
            "Position",
            "PFAM",
            "PANTHER",
            "KOG", 
            "GO",
            "Descr.",
            "AT def.",
            "OS def.", 
            "% similrty."
        };
        for(int i=0; i< data.length; i++) if (pieces.length > data[i]) {
            stb.append("<tr><td>")
                    .append(info[i])
                    .append("</td><td>")
                    .append(pieces[data[i]])
                    .append("</td></tr>");
        }
        stb.append("</table>");
        return stb.toString();
    }

    private static final int LINE_LENGTH = 80;
    public String toFasta() {
        StringBuilder ret = new StringBuilder(">"+name+"\t"+gff3EntryGene+"\n");
        String nucleotidos = gene;
        while(! nucleotidos.isEmpty()) {
            ret.append(nucleotidos.substring(0, Math.min(LINE_LENGTH, nucleotidos.length()))).append("\n");
            nucleotidos = nucleotidos.substring(Math.min(LINE_LENGTH, nucleotidos.length()));
        }
        return ret.toString();
    }

}
