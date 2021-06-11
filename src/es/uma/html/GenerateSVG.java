package es.uma.html;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/* This program generates a beautiful diagram with one column per Brachypodium variety.
 * Diagram is in SVG format, it is interactive and can be open in a browser.
 * Dehydrins are coloured depending of their relationship through a N/W similarity.
 * Genes may have different shapes:
 * - Square: It is marked as PF00257 and the chromosome is known.
 * - Triangle: It is not marked as PF00257 but contains the word "dehydrin" and the chromosome is known.
 * - Circle: It belongs to an unknown chromosome. Can be moved with the mouse when displayed in a browser.
 * - Diamond: It has been manually curated.
 */
public class GenerateSVG {
    private static final String INPUT_DATA_FILE = "datafiles/dehydrins_ecotypes.data"; // "G:\\Data_JIC\\NetBeansProjects\\Utilities\\src\\brachypodium\\varieties\\searchgenes\\dehydrins.sal";
    // Lengths of chromosomes taken from:
    // Genome sequencing and analysis of the model grass Brachypodium distachyon (https://www.nature.com/articles/nature08747)
    private static final int[] CHR_LENGTHS = {74834646, 59328898, 59892396, 48648102, 28648102};
    public static List<Variety> varieties = new ArrayList<>();
    private static String[] backgrounds = {};
    private static String[] titles = {};
    static {
        String line = null;
        try (BufferedReader in = Files.newBufferedReader(Paths.get("datafiles/Order.txt"));){
            List<String> colors = new ArrayList<>();
            List<String> texts = new ArrayList<>();
            while ((line = in.readLine()) != null)
                if (! line.startsWith("#")) {
                    colors.add(line.split("\t")[1]);
                    texts.add(line.split("\t")[2]);
                }
            backgrounds = colors.toArray(backgrounds);
            titles = texts.toArray(titles);
        } catch (Exception x) { System.err.println(line); System.err.println(Paths.get(".").toAbsolutePath().toString()); x.printStackTrace(); }
    }
    public static void main(String[] args) {
        loadGenes(); // This fills the variable: varieties.
        // verResumen();
        //Collections.sort(varieties);
        try (PrintStream myOut = new PrintStream(new BufferedOutputStream(new FileOutputStream("result.html")), true)) {
            myOut.println(initHTML());
            myOut.println(firstRow());
            FillColors.iterateVarietiesAgainstFirst(varieties); // This groups on Bdhn1 to Bdhn10
            // iterateVarietiesAgainstFirst(variedades); // This assigns reference Bdhns.
            for(int chrNum: new int[]{1, 2, 3, 4, 5}) {
                FillColors.iterateChromosomes(chrNum, varieties); // Model distachyon is at pos 0
                FillColors.iterateUnknownChromosomes(chrNum, varieties);
                myOut.println(initRow("chr"+chrNum));
                for(Variety v: varieties){
                    displayChr(myOut, v, chrNum);
                }
                myOut.println(endRow());
            }
            //
            myOut.println("<tr><td></td>");
            for(Variety v: varieties){
                displaySummary(myOut, v);
            }
            myOut.println("</tr>");
            //
            myOut.println("<tr><td></td>");
            /*  
            * Here we want to show the distribution of the lengths of the genes.
            * To normalize the distribution, we use as left and right limits, those of the 
            * model Brachypodium (minimum and maximum lengths, respectively). 
            * For any variety, any gene lower than the minimum or greater than the maximum
            * it will appear in a single left or right block.
            */
            List<Integer> sizes = Stream.of(varieties.get(0).genes).map(g -> g.protein.length()).collect(Collectors.toList());
            Collections.sort(sizes);
            int minimum = sizes.get(0);
            int maximum = sizes.get(sizes.size()-1) + 1;
            // To know the maximum size of a bar, we call the function as a stub
            int scale = 0;
            for(Variety v: varieties)
                scale = Integer.max(scale, displayDistribution(System.out, v, minimum, maximum, 100));
            // Once the left and right limits have been calculated.
            // And once the maximum height has been calculated (scale).
            // We can draw the graph
            for(Variety v: varieties)
                displayDistribution(myOut, v, minimum, maximum, Integer.min(18, scale));
            myOut.println("</tr>");
            //
            myOut.println(endHTML());
            System.out.println();
            showClustersOfGenes(varieties); // This call makes sense only when FillColors used do not group on Bdhn1 to Bdhn10
            //checkSegmentK(varieties);
            // Distribution.showDistributionBdhn1To10(varieties, titles); // This makes sense only when genes has been clasified over Bdhn1 to Bdhn10
        } catch (Exception x) { x.printStackTrace(); }
    }

    private static final int GENE_HEIGHT = 20;
    private static void displayChr(PrintStream myOut, Variety v, int chrNumber) {
        int chrLength = CHR_LENGTHS[chrNumber-1]/100_000;
        double scale = CHR_LENGTHS[chrNumber-1]/(double)v.chrLength[chrNumber-1];
        String keyword = "ALEA JACTA EST";
        //
        // REMOVE background-color
        v.bakcground = "rgb(255,255,255)";
        String text = "<td id="+chrNumber+" style='background-color:"+v.bakcground+";'><svg  onload=\"makeDraggable(evt)\" width="+(GENE_HEIGHT+4)+" height="+(chrLength+GENE_HEIGHT+10+20)+">"
                //+ "<rect x=5 y=5 width=10 height="+(chrLength+GENE_HEIGHT)+" style='fill:rgb(248,222,166);stroke-width:2;stroke:rgb(0,0,0)' />"
                //+ "<text x=750 y=13 fill='black' text-anchor='middle'>"+init+"</text>"
                //+ "<circle cx="+(chrLength/2)+" cy=20 r=6 stroke='black' stroke-width=2 fill='red' />"
                + keyword
                + "</svg></td>";
        String genesHTML = "";
        int lastY = -1;
        for(Gene g: v.retrieveFilteredGenesinOrder(chrNumber)) {
            int posY = (int)(scale*g.init/100_000+5);
            if (posY < lastY) posY = lastY;
            lastY = posY + GENE_HEIGHT;
            String points;
            switch (g.shape) {
                case "RECT":
                    genesHTML += "<rect id=\""+g.getHTML()+"\" x=2 y="+(posY)+" width="+GENE_HEIGHT+" height="+GENE_HEIGHT+" style='fill:"+g.getColor()+";stroke-width:2;stroke:rgb(0,0,0)' onmousemove='showTooltip(evt);' onmouseout='hideTooltip(evt);' ></rect>";
                    break;
                case "CIRCLE":
                    genesHTML += "<circle class=\"draggable\" id=\""+g.getHTML()+"\" cx="+(2+GENE_HEIGHT/2)+" cy="+(posY+GENE_HEIGHT/2)+" r="+(GENE_HEIGHT/2)+" style='fill:"+g.getColor()+";stroke-width:2;stroke:rgb(0,0,0)' onmousemove='showTooltip(evt);' onmouseout='hideTooltip(evt);' ></circle>";
                    break;
                case "DIAMOND":
                    points = "12,"+posY+    ","+(12+GENE_HEIGHT/2)+","+(posY+GENE_HEIGHT/2)+    ",12,"+(posY+GENE_HEIGHT)+   ","+(12-GENE_HEIGHT/2)+","+(posY+GENE_HEIGHT/2);
                    genesHTML += "<polygon  id=\""+g.getHTML()+"\" points="+points+" style='fill:"+g.getColor()+";stroke-width:2;stroke:rgb(0,0,0)' onmousemove='showTooltip(evt);' onmouseout='hideTooltip(evt);' ></polygon >";
                    break;
                default:
                    points = "12,"+posY+","+(12+GENE_HEIGHT/2)+","+(posY+GENE_HEIGHT)+","+(12-GENE_HEIGHT/2)+","+(posY+GENE_HEIGHT);
                    genesHTML += "<polygon  id=\""+g.getHTML()+"\" points="+points+" style='fill:"+g.getColor()+";stroke-width:2;stroke:rgb(0,0,0)' onmousemove='showTooltip(evt);' onmouseout='hideTooltip(evt);' ></polygon >";
                    break;
            }
        }
        myOut.println(text.replace(keyword, genesHTML));
    }
    
    public static void loadGenes() {
        try (BufferedReader in = Files.newBufferedReader(Paths.get(INPUT_DATA_FILE));){
            int backgroundPos = 0;
            while (loadVariety(in, backgroundPos++));
        } catch (Exception x) { x.printStackTrace(); }
    }

    private static boolean loadVariety(BufferedReader in, int backgroundPos) throws IOException {
        String name = in.readLine();
        // System.out.println(name);
        if (name == null) return false;
        name = name.replace(".annotation_info.txt", "");
        int longitud = Integer.parseInt(in.readLine().split("\t")[1]);
        Variety v = new Variety(name, longitud, (backgroundPos<backgrounds.length)?backgrounds[backgroundPos]:"rgb(255,255,255)");
        for(int i=0; i<5; i++)
            v.chrLength[i] = Integer.parseInt(in.readLine().split("\t")[1]);
        for(int i=0; i< longitud; i++) {
            v.genes[i] =  new Gene(in.readLine());
            v.genes[i].protein = in.readLine();
            v.genes[i].gff3EntryGene = in.readLine();
            v.genes[i].gene = in.readLine();
        }
        int current = v.genes.length;
        v.removeGenesWithoutFasta();
        if (current != v.genes.length) System.out.println("Removed: "+ (current - v.genes.length));
        current = v.genes.length;
        if (current != v.genes.length) System.out.println("Removed: "+ (current - v.genes.length));
        varieties.add(v);
        return true;
    }
    
         private static final String TABLE_BDHN = 
                "<td style=\"border:0px;\">&nbsp;</td><td style=\"border:0px;vertical-align:top\">\n" +
                "<table border=1 style='font-size:18px;font-family:Verdana'>\n" +
                "<tr><th>DHN Colour</th><th>B. distachyon DHN</th><th>Bd1</th><th>Bd2</th><th>Bd3</th><th>Bd4</th><th>Bd5</th></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(102,204,255);stroke-width:0'></svg></td><td align=\"center\">Bdhn1</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>2</td><td>51</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(204,204,255);stroke-width:0'></svg></td><td align=\"center\">Bdhn2</td><td>&nbsp;</td><td>&nbsp;</td><td>20</td><td>34</td><td>&nbsp;</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(255,204,255);stroke-width:0'></svg></td><td align=\"center\">Bdhn3</td><td>54</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(255,204,204);stroke-width:0'></svg></td><td align=\"center\">Bdhn4</td><td>&nbsp;</td><td>&nbsp;</td><td>27</td><td>19</td><td>&nbsp;</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(255,204,153);stroke-width:0'></svg></td><td align=\"center\">Bdhn5</td><td>&nbsp;</td><td>&nbsp;</td><td>34</td><td>19</td><td>&nbsp;</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(255,255,153);stroke-width:0'></svg></td><td align=\"center\">Bdhn6</td><td>&nbsp;</td><td>&nbsp;</td><td>6</td><td>44</td><td>2</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(204,255,153);stroke-width:0'></svg></td><td align=\"center\">Bdhn7</td><td>&nbsp;</td><td>&nbsp;</td><td>24</td><td>30</td><td>&nbsp;</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(153,255,153);stroke-width:0'></svg></td><td align=\"center\">Bdhn8</td><td>&nbsp;</td><td>&nbsp;</td><td>22</td><td>33</td><td>&nbsp;</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(154,254,214);stroke-width:0'></svg></td><td align=\"center\">Bdhn9</td><td>&nbsp;</td><td>49</td><td>&nbsp;</td><td>&nbsp;</td><td>10</td></tr>\n" +
                "<tr style='line-height:12px;'><td style='padding:0;margin:0;'><svg width=120 height=24><rect id='Bdhn1' x=0 y=0 width=190 height=24 style='fill:rgb(102,255,255);stroke-width:0'></svg></td><td align=\"center\">Bdhn10</td><td>52</td><td>&nbsp;</td><td>&nbsp;</td><td>2</td><td>&nbsp;</td></tr>\n" +
                "</table>\n" +
                "</td>\n" +
                "</td>";
        
        private static final String HTMLHeader = "<!doctype html>\n" +
        "<html lang=\"en\">\n" +
        "<head>\n" +
        "  <meta charset=\"utf-8\">\n" +
        "  <title>DEHYDRINS HTML5 results</title>\n" +
        "  <meta name=\"description\" content=\"HTML5 results\">\n" +
        "  <meta name=\"author\" content=\"SitePoint\">\n" +
            "<style>"+
                "#tooltip {\n" +
                    "  background: cornsilk;\n" +
                    "  border: 1px solid black;\n" +
                    "  border-radius: 5px;\n" +
                    "  padding: 5px;\n" +
                    "}"+
                "table {\n" +
                "  border-collapse: collapse;\n" +
                "}\n" +
                "table td, table th {\n" +
                "  border: 1px solid black;\n" +
                "  border-right: 1px solid #ddd;\n" +
                "  border-left: 1px solid #ddd;\n" +
                "}\n" +
                "table tr:first-child th {\n" +
                "  border-top: 1px solid black;\n" +
                "}\n" +
                "table tr td:first-child,\n" +
                "table tr th:first-child {\n" +
                "  border-left: 1px solid black;\n" +
                "  border-right: 1px solid black;\n" +
                "}\n" +
                "table tr td:last-child,\n" +
                "table tr th:last-child {\n" +
                "  border-right: 1px solid black;\n" +
                "}\n" +
                ".draggable {\n" +
                "  cursor: move;\n" +
                "}" +
            "</style><script>"+
                "function showTooltip(evt, text) {\n" +
            "  let tooltip = document.getElementById(\"tooltip\");\n" +
            "  tooltip.innerHTML = evt.target.id;\n" +
            "  tooltip.style.display = \"block\";\n" +
            "  tooltip.style.left = evt.pageX + 10 + 'px';\n" +
            "  tooltip.style.top = evt.pageY + 10 + 'px';\n" +
            "  let celdas = evt.target.parentNode.parentNode.parentNode.getElementsByTagName(\"TD\");\n" +
            "  for(var i=1;i<celdas.length;i++) {\n" +
            "		let rects = celdas[i].childNodes[0].childNodes;\n" +
            "		for(var j=0;j<celdas.length;j++){\n" +
            "			if (rects[j] == null) continue;\n" +
            "			if (rects[j].style.fill == evt.target.style.fill)\n" +
            "				rects[j].style.stroke = \"rgb(255,0,0)\";\n" +
            "		}\n" +
            "  }"+
            "}\n" +
            "\n" +
            "function hideTooltip(evt) {\n" +
            "  var tooltip = document.getElementById(\"tooltip\");\n" +
            "  tooltip.style.display = \"none\";\n" +
            "  let celdas = evt.target.parentNode.parentNode.parentNode.getElementsByTagName(\"TD\");\n" +
            "  for(var i=1;i<celdas.length;i++) {\n" +
            "		let rects = celdas[i].childNodes[0].childNodes;\n" +
            "		for(var j=0;j<celdas.length;j++){\n" +
            "			if (rects[j] == null) continue;\n" +
            "			if (rects[j].style.fill == evt.target.style.fill)\n" +
            "				rects[j].style.stroke = \"rgb(0,0,0)\";\n" +
            "		}\n" +
            "  }"+
            "}\n"+
            "var selectedElement = false;\n" +
            "function makeDraggable(evt) {\n" +
            "  var svg = evt.target;\n" +
            "  svg.addEventListener('mousedown', startDrag);\n" +
            "  svg.addEventListener('mousemove', drag);\n" +
            "  svg.addEventListener('mouseup', endDrag);\n" +
            "  svg.addEventListener('mouseleave', endDrag);\n" +
            "  function getMousePosition(evt) {\n" +
            "	  var CTM = svg.getScreenCTM();\n" +
            "	  return {\n" +
            "		x: (evt.clientX - CTM.e) / CTM.a,\n" +
            "		y: (evt.clientY - CTM.f) / CTM.d\n" +
            "	  };\n" +
            "  }\n" +
            "  function startDrag(evt) {\n" +
            "	  if (evt.target.classList.contains('draggable')) {\n" +
            "		selectedElement = evt.target;\n" +
            "	  }\n" +
            "  }\n" +
            "  function drag(evt) {\n" +
            "	  if (selectedElement) {\n" +
            "		evt.preventDefault();\n" +
            "		var coord = getMousePosition(evt);\n" +
            "		selectedElement.setAttributeNS(null, \"cy\", coord.y);\n" +
            "	  }\n" +
            "  }\n" +
            "  function endDrag(evt) {\n" +
            "	selectedElement = null;\n" +
            "  }\n" +
            "}"+
            "</script>"+
        "</head>\n" +
        "<body>"+
        "<h1 align='center' >DEHYDRINS in 54 ecotypes of Brachypodium distachyon<br/></h1>"+
        "<div id='tooltip' display='none' style='position: absolute; display: none;'></div>"+
        "<table style=\"border:0px;border-collapse:separate;\"><tr style=\"border:0px;\"><td style=\"border:0px;\">\n" +
//        "<table>\n" +
//        "<tr><th>Colour</th><th>Meaning</th></tr>\n" +
//        "<tr><td style='background-color:rgb(220,255,255);'>&nbsp;</td><td>Main Reference</td></tr>\n" +
//        "<tr><td style='background-color:rgb(255,220,255);'>&nbsp;</td><td>Tolerant</td></tr>\n" +
//        "<tr><td style='background-color:rgb(255,240,220);'>&nbsp;</td><td>Intermediary</td></tr>\n" +
//        "<tr><td style='background-color:rgb(255,255,220);'>&nbsp;</td><td>Susceptible</td></tr>\n" +
//        "<tr><td style='background-color:rgb(204,255,204);'>&nbsp;</td><td>PWC&gt;70%</td></tr>\n" +
//        "<tr><td style='background-color:rgb(230,255,230);'>&nbsp;</td><td>50%&gt;PWC&gt;70%</td></tr>\n" +
//        "<tr><td style='background-color:rgb(255,255,153);'>&nbsp;</td><td>PWC&lt;55%</td></tr>\n" +
//        "<tr><td style='background-color:rgb(204,255,255);'>&nbsp;</td><td>Dry medierranean lands</td></tr>\n" +
//        "<tr><td style='background-color:rgb(204,255,229);'>&nbsp;</td><td>Olive and cereal lands</td></tr>\n" +
//        "<tr><td style='background-color:rgb(229,255,204);'>&nbsp;</td><td>Subhumid forests and meadows</td></tr>\n" +
//        "</table>\n" +
        "</td><td style=\"border:0px;\">&nbsp;</td><td style=\"border:0px;vertical-align:top\">\n" +
        "<table style='font-size:18px;font-family:Verdana'>\n" +
        "<tr><th>Shape</th><th>Meaning</th></tr>\n" +
        "<tr><td align='center'><svg width='24px' height='24px'><rect x=2 y=2 width=20 height=20 style='fill:#FFFFFF;stroke-width:2;stroke:rgb(0,0,0)'></rect></svg></td><td>Annotated as PF0257</td></tr>\n" +
        "<tr><td align='center'><svg width='24px' height='24px'><polygon points=12,2,2,22,22,22 style='fill:#FFFFFF;stroke-width:2;stroke:rgb(0,0,0)'></polygon></svg></td><td>Only similar to other specie's dehydrin</td></tr>\n" +
        "<tr><td align='center'><svg width='24px' height='24px'><polygon points=12,2,2,12,12,22,22,12 style='fill:#FFFFFF;stroke-width:2;stroke:rgb(0,0,0)'></polygon></svg></td><td>Manually found</td></tr>\n" +
        "<tr><td align='center'><svg width='24px' height='24px'><circle cx=12 cy=12 r=10 style='fill:#FFFFFF;stroke-width:2;stroke:rgb(0,0,0)'></rect></svg></td><td>Belongs to an unknown chromosome</td></tr>\n" +
        "</table>\n" +
        "</td>"+
        TABLE_BDHN +
        "</tr></table>\n" +
        "<p>&nbsp;</p>"+
        "<table border=1>";
    public static String initHTML() {
        return HTMLHeader;
    }
    
    private static final String HTMLFooter = "</table><p>&nbsp;</p><p>&nbsp;</p><p>&nbsp;</p><p>&nbsp;</p></body>\n" +
    "</html>";
    private static String endHTML() {
        return HTMLFooter;
    }

    private static String initRow(String nameRow) {
        return "<tr><td style='text-anchor:end;font-family:Verdana;font-size:20px;'><strong>"+nameRow+"</strong></td>";
    }
    private static String endRow() {
        return "</tr>";
    }

    private static String firstRow() {
        StringBuilder stb = new StringBuilder();
        stb.append("<tr><th></th>");
        for(String title: titles) {
            stb.append("<th><svg width='24px' height='120px'>")
                    .append("<text transform='rotate(270)' style='text-anchor:end;font-family:Verdana;font-size:20px;' y='19'>")
                    .append(title)
                    .append("</text></svg></th>");
        }
        stb.append("</tr>");
        return stb.toString();
    }

    private static void displaySummary(PrintStream myOut, Variety v) {
        int unk = 0;
        for(Gene g: v.genes){
            if (g.getChrNumber() == -1) unk++;
            //if (belongToUnknown(g)) unk++;
        }
        myOut.print("<td><svg width='24px' height='70px'><text transform='rotate(270)' style='text-anchor:end;font-family:Verdana;font-size:20;' style='text-anchor:end;' y='19'>"+unk+"("+v.genes.length+")</text></svg></td>");
    }
    
    private static int displayDistribution(PrintStream myOut, Variety v, int first, int last, int scale) {
        int NUM_BARS = 4;
        List<Integer> sizes = Stream.of(v.genes).map(g -> g.protein.length()).collect(Collectors.toList());
        myOut.print("<td><svg width='20px' height='92px'>"); // <text transform='rotate(270)' style='text-anchor:end;' y='15'>");
        int numGenes = (int) sizes.stream().filter(x -> (x<first)).count();
        int maxNumGenes = numGenes;
        int rectSize = numGenes*18/scale;
        myOut.print("<rect id="+numGenes+" x="+(19-rectSize)+" y=76 width="+rectSize+" height=15 style='fill:darkBlue;stroke-width:0;stroke:rgb(0,0,0)' ></rect>");
        for(int i=0; i<NUM_BARS; i++){
            int start = first + i * (last - first) / NUM_BARS;
            int end = first + (i+1) * (last - first) / NUM_BARS;
            numGenes = (int) sizes.stream().filter(x -> (x>=start) && (x<end)).count();
            maxNumGenes = Integer.max(maxNumGenes, numGenes);
            rectSize = numGenes*18/scale;
            myOut.print("<rect id="+numGenes+" x="+(19-rectSize)+" y="+(76-15*(i+1))+" width="+rectSize+" height=15 style='fill:darkBlue;stroke-width:0;stroke:rgb(0,0,0)' ></rect>");
        }
        numGenes = (int) sizes.stream().filter(x -> (x>=last)).count();
        maxNumGenes = Integer.max(maxNumGenes, numGenes);
        rectSize = numGenes*18/scale;
        myOut.print("<rect id="+numGenes+" x="+(19-rectSize)+" y="+(76-15*(NUM_BARS+1))+" width="+rectSize+" height=15 style='fill:darkBlue;stroke-width:0;stroke:rgb(0,0,0)' ></rect>");
        myOut.print("</svg></td>");
        return maxNumGenes;
    }
    //
    // This piece of code is very similar to that of the same name in FillColors.
    //
    static void iterateVarietiesAgainstFirst(List<Variety> variedades) {
        final int[] POSITIONS = {10, 3, 9, 8, 7, 2, 6, 4, 5, 1};
        int index = -1;
        // A different colour is assigned to each gene of the model (first variety, usually B. distachyon).
        Variety vModelo = variedades.get(0);
        for(Gene g: vModelo.genes){
            g.dhnCode = ""+POSITIONS[++index];
        }
        // For any other gene of any other variety, we assign to it the reference
        // of the most similar gene in the model.
        for(int i=1; i<variedades.size(); i++) {
            for(Gene g: variedades.get(i).genes){
                Gene geneModelo = FillColors.mostSimilarGene(g, vModelo);
                g.dhnCode = geneModelo.dhnCode;
            }
        }
    }

    private static final String DIR_CLUSTERS = "H:\\triticum\\Brachypodium\\ClusterIntraespecie\\";
    private static final StringBuilder sbResto = new StringBuilder();
    private static void showClustersOfGenes(List<Variety> variedades) {
        int clusterNum = 1;
        //System.out.println("*************************************************************************");
        //System.out.println("*************************************************************************");
        List<Gene> fase2 = new ArrayList<>();
        for(int chrNum: new int[]{1, 2, 3, 4, 5}) {
            for (Variety v: variedades) {
                for (int i=0; i<v.genes.length; i++) if (v.genes[i] != null && v.genes[i].getChrNumber() == chrNum) {
                    String filename = DIR_CLUSTERS+"Cluster_"+String.format("%03d", clusterNum++)+".fasta";
                    showAndRemoveGenes(variedades, v.genes[i], chrNum, filename);
                    //g.color = null;
                }
            }
        }
        // Display genes of other pseudomolecules.
        for (Variety v: variedades)
            System.out.print("\t"+v.name);
        System.out.println("\nUnclustered");
        for (Variety v: variedades) {
            System.out.print("\t");
            for (int i=0; i<v.genes.length; i++) if (v.genes[i] != null) {
                fase2.add(v.genes[i]);
                System.out.printf("%s(%d)",v.genes[i].dhnCode, v.genes[i].getChrNumber());
                v.genes[i] = null;
            }
        }
        int index = 0;
        clusterNum = 100;
        System.out.println("\nGenes from other pseudomolecules "+fase2.size());
        for (Gene g: fase2) if (g.color == null) {
            g.color = FillColors.COLORS[index];
            System.out.println(g.toFasta());
            double optimalScore = (double)FillColors.needlemanWunsch(g.protein, g).score;
            int counter = 1;
            StringBuilder sb = new StringBuilder();
            sb.append(g.toFasta()).append("\n");
            for (Gene g2: fase2) if (g2.color == null && !g.protein.equals("")) {
                FillColors.NWResult result = FillColors.needlemanWunsch(g.protein, g2);
                if (result.score/optimalScore > 0.95) {
                    result.gene.color = FillColors.COLORS[index];
                    System.out.println(g2.toFasta());
                    sb.append(g2.toFasta()).append("\n");
                    counter++;
                }
            }
            if (counter>2) try { Files.write(Paths.get(DIR_CLUSTERS+"Cluster_"+String.format("%03d", clusterNum++)+".fasta"), sb.toString().getBytes()); } catch (Exception x) { x.printStackTrace(); }
            else sbResto.append(sb);
            System.out.println(counter + " *************************************************************************");
            index++;
        }
        try { Files.write(Paths.get(DIR_CLUSTERS+"Resto.fasta"), sbResto.toString().getBytes()); } catch (Exception x) { x.printStackTrace(); }
    }

    private static void showAndRemoveGenes(List<Variety> variedades, Gene gene, int chrNum, String filename) {
        int counter = 0;
        StringBuilder sb = new StringBuilder();
        for (Variety v: variedades)
            System.out.print("\t"+v.name);
        System.out.println();
        System.out.print(filename);
        for (Variety v: variedades) {
            System.out.print("\t");
            for (int i=0; i<v.genes.length; i++) if (v.genes[i] != null && v.genes[i].getChrNumber() == chrNum && gene.color.equals(v.genes[i].color) ) {
                //System.out.println(v.genes[i].toFasta());
                sb.append(v.genes[i].toFasta()).append("\n");
                counter++;
                System.out.printf("%s(%d)",v.genes[i].dhnCode, v.genes[i].getChrNumber());
                v.genes[i] = null;
            }
        }
        System.out.println();
        if (counter>2) try { Files.write(Paths.get(filename), sb.toString().getBytes()); } catch (Exception x) { x.printStackTrace(); }
        else sbResto.append(sb);
        //System.out.println("*************************************************************************");
    }

}

