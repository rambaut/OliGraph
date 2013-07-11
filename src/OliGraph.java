import jebl.evolution.align.*;
import jebl.evolution.align.scores.*;

import java.io.*;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class Critters {

    public Critters(String[] args) throws IOException {

        int minLength = Integer.parseInt(args[0]);

        String fileName1 = args[1];
        Map<String, String> readMap1 = readFasta(fileName1, minLength);

        String fileName2;
        Map<String, String> readMap2;


        final PrintWriter outLog = new PrintWriter(new FileWriter("alignments.txt"));
        final PrintWriter outEdges = new PrintWriter(new FileWriter("edges.txt"));
        outEdges.println("from\tscore\tto\titeration\tmatches\tmismatches\tgaps\tidentity\tdistance");

        final PrintWriter outNodes = new PrintWriter(new FileWriter("nodes.txt"));
        outNodes.println("name\ttime\tcount\tlineages1\tlineages2\tlineages3\tdescendants1\tdescendants2\tdescendants3");

        Scores scores = new NucleotideScores(5,-4);

        float d = 20.0f;
        float e = 1.0f;

//        Align align = new NeedlemanWunschLinearSpaceAffine(scores, d, e);
//        align = new OldNeedlemanWunschAffine(scores, d, e);
        Align align = new SmithWatermanLinearSpaceAffine(scores, d, e);

        List<String> names = new ArrayList<String>();
        names.addAll(readMap1.keySet());

        Map<String, Set<String>> descendants = new HashMap<String, Set<String>>();

        for (int i = 2; i < args.length; i++) {
            fileName2 = args[i];
            readMap2 = readFasta(fileName2, minLength);

            names.addAll(readMap2.keySet());

            outLog.println("File: " + fileName1 + " vs File: " + fileName2);
            for (String name2 : readMap2.keySet()) {
                String sequence2 = readMap2.get(name2);
                String id2 = getID(name2);

                float bestScore = 0;
                String bestKey = null;

                for (String name1 : readMap1.keySet()) {
                    String sequence1 = readMap1.get(name1);

                    align.doAlignment(sequence1, sequence2);
                    float score = align.getScore();
                    if (score > bestScore) {
                        bestScore = score;
                        bestKey = name1;
                    }
                }

                Set<String> links = descendants.get(bestKey);
                if (links == null) {
                    links = new HashSet<String>();
                    descendants.put(bestKey, links);
                }
                links.add(name2);

                String sequence1 = readMap1.get(bestKey);
                align.doMatch(new Output() {
                    public final void print(String s) {
                        outLog.print(s);
                    }

                    public final void println(String s) {
                        outLog.println(s);
                    }

                    public final void println() {
                        outLog.println();
                    }
                }, name2 + " -> " + bestKey);
                align.doAlignment(sequence1, sequence2);

                String[] match = align.getMatch();
                int[] counts = matchCounts(match);

                int matches = counts[0];
                int mismatches = counts[1];
                int gaps = counts[2];
                double identity = (double)counts[0] / match[0].length();

                String bestId = getID(bestKey);

                outEdges.println(id2 + "\t" + align.getScore() + "\t" + bestId + "\t" + i +
                        "\t" + matches + "\t" + mismatches + "\t" + gaps + "\t" + identity );
            }

            System.out.println();
            System.out.println();
            System.out.println();

            readMap1 = readMap2;
            fileName1 = fileName2;
        }

        for (String name : names) {
            int d1 = 0;
            int d2 = 0;
            int d3 = 0;

            int c1 = 0;
            int c2 = 0;
            int c3 = 0;

            String id = getID(name);
            int time = getTime(name);
            int count = getCount(name);

            Set<String> descendants1 = descendants.get(name);
            if (descendants1 != null) {
                d1 = descendants1.size();
                for (String name2 : descendants1) {
                    c1 += getCount(name2);

                    Set<String> descendants2 = descendants.get(name2);
                    if (descendants2 != null) {
                        d2 += descendants2.size();
                        for (String name3 : descendants2) {
                            c2 += getCount(name3);

                            Set<String> descendants3 = descendants.get(name3);
                            if (descendants3 != null) {
                                d3 += descendants3.size();
                                for (String name4 : descendants3) {
                                    c3 += getCount(name4);
                                }
                            }

                        }
                    }
                }
            }
            outNodes.println(id + "\t" + time + "\t" + count + "\t" + d1 + "\t" + d2  + "\t" + d3 + "\t" + c1 + "\t" + c2  + "\t" + c3 );
        }

        outLog.close();
        outEdges.close();
        outNodes.close();
    }

    private String getID(String name) {
        String[] parts = name.split("_");
        return parts[3] + "_" + parts[0];
    }

    private int getTime(String name) {
        String[] parts = name.split("_");
        return Integer.parseInt(parts[3]);
    }

    private int getCount(String name) {
        String[] parts = name.split("_");
        return Integer.parseInt(parts[5]);
    }

    int[] matchCounts(String[] match) {
        int[] matchCounts = new int[3];
        for (int i = 0; i < match[0].length(); i++) {
            char c1 = match[0].charAt(i);
            char c2 = match[1].charAt(i);

            if (c1 == c2) {
                matchCounts[0] += 1;
            } else if (c1 != '-' && c2 != '-') {
                matchCounts[1] += 1;
            } else {
                matchCounts[2] += 1;
            }
        }
        return matchCounts;
    }

    public Map<String, String> readFasta(String inFileName, int minLength) {
        int count = 0;
        Map<String, String> sequences = new HashMap<String, String>();

        try {
            BufferedReader reader = new BufferedReader(new FileReader(inFileName));

            String line = reader.readLine();
            while (line != null) {
                if (line.startsWith(">")) {
                    String label = line.substring(1).split(" ")[0];

                    StringBuilder seq = new StringBuilder();
                    line = reader.readLine();
                    if (line == null || line.startsWith(">")) {
                        throw new RuntimeException("Sequence missing.");
                    }

                    while (line != null && !line.startsWith(">")) {
                        for (char c : line.toCharArray()) {
                            if (!Character.isWhitespace(c)) {
                                seq.append(c);
                            }
                        }
                        line = reader.readLine();
                    }

                    String seqString = seq.toString();
                    if (seqString.length() >= minLength) {
                        sequences.put(label, seqString);
                    } else {
                        count ++;
                    }
                } else {
                    line = reader.readLine();
                }
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return null;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }

        System.err.println("Read file, " + inFileName + ": " + sequences.keySet().size() + " sequences (" + count + " with length < " + minLength + " skipped)");
        return sequences;
    }

    public static void main(String[] args) {

        if (args.length >= 3) {
            try {
                new Critters(args);
            } catch (IOException e) {
                e.printStackTrace();
            }
        } else {
            System.err.println("Usage: critters <min_length> <infilename1> <infilename2> [<infilename3 ...]");
            return;
        }
    }
}
