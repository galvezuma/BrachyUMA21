package es.uma.motif;

/**
* This class allows implementing a motif searcher using an exact pattern. It allows a
* particular number of mismatches. Usually, the number of allowed mismatches depends on
* the actual length of the exact pattern.
*  
* @author galvez
*/
public class PerfectMatchArgument {
   protected String aminoAcidSequence;
   protected String name;
   protected int numMismatches;
   public PerfectMatchArgument(String aminoAcidSequence, String name, int numMismatches) {
       this.aminoAcidSequence = aminoAcidSequence;
       this.name = name;
       this.numMismatches = numMismatches;
   }
}
