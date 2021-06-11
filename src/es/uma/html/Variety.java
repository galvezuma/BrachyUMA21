package es.uma.html;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Predicate;

public class Variety implements Comparable<Variety> {
    public String name;
    int[] chrLength;
    public Gene[] genes;
    String bakcground;

    Variety(String name, int longitud, String background) {
        this.name = name;
        this.genes = new Gene[longitud];
        this.chrLength = new int[5];
        this.bakcground = background;
    }
    
    public void removeGenesWithoutFasta(){
        filterGenes(g -> g.protein != null && ! g.protein.equals(""));
    }
    
    public void filterGenes(Predicate<Gene> filtro) {
        List<Gene> list = new LinkedList<>();
        Arrays.asList(genes).stream().filter(filtro).forEach(list::add);
        genes = new Gene[list.size()];
        list.toArray(genes);
    }
    
    @Override
    public int compareTo(Variety v2) {
            return this.genes.length - v2.genes.length;
    }

    Iterable<Gene> retrieveFilteredGenesinOrder(int chrNumber) {
        List<Gene> ret = new ArrayList<>();
        for(Gene g: genes) if (g.getChrNumber() == chrNumber)
            ret.add(g);
        ret.sort(null);
        return ret;
    }
}
