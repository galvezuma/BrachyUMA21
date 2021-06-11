package es.uma.motif;

import java.util.regex.Pattern;


/**
* This class allows implementing a motif searcher using a regular expression.
*  
* @author galvez
*/
public class PatternArgument {
    protected Pattern pattern;
    protected String name;
    public PatternArgument(Pattern pattern, String name) {
        this.pattern = pattern;
        this.name = name;
    }
}
