import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class NeedlemanWunsch {

    static Map<String, String> map = new HashMap<>();

    static File sequenceOneFile;
    static File sequenceTwoFile;
    static File output;

    static int matchScore = 1;
    static int mismatchPenalty = -1;
    static int regularGapPenalty = -2;
    static int initialGapPenalty = regularGapPenalty;

    static int matchCount = 0;
    static int misMatchCount = 0;
    static int gapCount = 0;

    static boolean convertToAA = false;
    static boolean trimToStartStop = false;

    public static void main(String[] args) throws Exception {
        if(args.length < 3) {
            throw new IllegalArgumentException("Received insufficient number of arguments. Expected usage: \"java NeedlemanWunsch sequence1 sequence2 output\"");
        }
        sequenceOneFile = new File("../data/contigs/" + args[0]);
        sequenceTwoFile = new File("../data/" + args[1]);
        output = new File("../data/alignments/" + args[2]);

        Scanner terminal = new Scanner(System.in);
        String input;

        //  Input files checks
        if(!sequenceOneFile.exists()) {
            System.out.println("File \"" + sequenceOneFile + "\" not found in data/contigs folder. Rerun program.");
            return;
        }
        if(!sequenceTwoFile.exists()) {
            System.out.println("File \"" + sequenceTwoFile + "\" not found in data/contigs folder. Rerun program.");
            return;
        }

//        //  Output file check
//        if(output.exists() && !output.isDirectory()) {
//            System.out.print(output + " already exists. Would you like to overwrite it? [y/n] ");
//            input = terminal.nextLine().toLowerCase();
//            while(!input.matches("^y|(yes)|n|(no)$")) {
//                System.out.print("Expected [y/n]: ");
//                input = terminal.nextLine().toLowerCase();
//            }
//            if(input.matches("^n|(no)$")) {
//                System.out.println("Terminating program. Select a new output location.");
//                return;
//            }
//        }

        //  Score and Penalty Handling
        if(args.length == 3) {
//            System.out.println("Received no match or penalty arguments. Proceeding with defaults:\n" +
//                    "\tMatch Score: " + matchScore + "\n" +
//                    "\tMismatch Penalty: " + mismatchPenalty + "\n" +
//                    "\tGap Penalty: " + regularGapPenalty + "\n" +
//                    "\tInitial Gap Penalty: " + initialGapPenalty + "\n");
        } else if(args.length != 6 && args.length != 7) {
//            System.out.println("Received invalid amount of arguments. Expected usages:\n" +
//                    "\"java NeedlemanWunsch sequence1 sequence2 output\" or\n" +
//                    "\"java NeedlemanWunsch sequence1 sequence2 output [Match Score] [Mismatch Penalty] [Gap Penalty]\" or\n" +
//                    "\"java NeedlemanWunsch sequence1 sequence2 output [Match Score] [Mismatch Penalty] [Gap Penalty] [Initial Gap Penalty]\"\n");
//            System.out.println("Terminating program.");
            return;
        } else {
            matchScore = Integer.parseInt(args[3]);
            mismatchPenalty = Integer.parseInt(args[4]);
            regularGapPenalty = Integer.parseInt(args[5]);
            if(args.length == 7) {
                initialGapPenalty = Integer.parseInt(args[6]);
            } else {
                initialGapPenalty = regularGapPenalty;
            }
//            System.out.println("Proceeding the following values:\n" +
//                    "\tMatch Score: " + matchScore + "\n" +
//                    "\tMismatch Penalty: " + mismatchPenalty + "\n" +
//                    "\tGap Penalty: " + regularGapPenalty + "\n" +
//                    "\tInitial Gap Penalty: " + initialGapPenalty + "\n");
        }

        initHashmap();

        //  NeedlemanWunsch Algorithm

        StringBuilder labelOne = new StringBuilder();
        StringBuilder sequenceOne = new StringBuilder();
        StringBuilder labelTwo = new StringBuilder();
        StringBuilder sequenceTwo = new StringBuilder();
        parseDnaFromFile(sequenceOneFile, labelOne, sequenceOne);
        parseDnaFromFile(sequenceTwoFile, labelTwo, sequenceTwo);

        if(convertToAA) {
            sequenceOne = dnaToAa(sequenceOne);
            sequenceTwo = dnaToAa(sequenceTwo);
        } else if(trimToStartStop) {
            sequenceOne = trimToStartAndStop(sequenceOne);
            sequenceTwo = trimToStartAndStop(sequenceTwo);
        }

        int lengthOne = sequenceOne.length();
        int lengthTwo = sequenceTwo.length();
        int[][] scoreMatrix = new int[lengthTwo + 1][lengthOne + 1];

        int i, j, match, sequenceOneGapScore, sequenceTwoGapScore;

        //  Alignment score matrix creation

        for(i = 0; i <= lengthTwo; i++) {  //  Populating first column
            scoreMatrix[i][0] = initialGapPenalty * i;
        }
        for(j = 0; j <= lengthOne; j++) {  //  Populating first row
            scoreMatrix[0][j] = initialGapPenalty * j;
        }

        for(i = 1; i <= lengthTwo; i++) {
            for(j = 1; j <= lengthOne; j++) {

                //  match = upper left diagonal + the match ( or mismatch ) score of the two characters
                match = scoreMatrix[i - 1][j - 1] + score(sequenceOne.charAt(j - 1), sequenceTwo.charAt(i - 1));

                //  Calculate gap scores
                sequenceOneGapScore = scoreMatrix[i - 1][j] + ( j == lengthOne ? initialGapPenalty : regularGapPenalty);
                sequenceTwoGapScore = scoreMatrix[i][j - 1] + ( i == lengthTwo ? initialGapPenalty : regularGapPenalty);

                //  Find greatest value, assign it to scoreMatrix[i][j]
                if(match >= Math.max(sequenceTwoGapScore, sequenceOneGapScore)) {
                    scoreMatrix[i][j] = match;
                } else if(sequenceOneGapScore >= Math.max(match, sequenceTwoGapScore)) {
                    scoreMatrix[i][j] = sequenceOneGapScore;
                } else {
                    scoreMatrix[i][j] = sequenceTwoGapScore;
                }
            }
        }
        StringBuilder modifiedSequenceOne = new StringBuilder();
        StringBuilder modifiedSequenceTwo = new StringBuilder();
        StringBuilder comparator = new StringBuilder();

        i = lengthTwo;
        j = lengthOne;

        //  Finding most optimal path
        while(i > 0 || j > 0) {

            //  If upper left ( match / mismatch ) is parent...
            if(i > 0 && j > 0 && scoreMatrix[i][j] == (scoreMatrix[i - 1][j - 1] + score(sequenceOne.charAt(j - 1), sequenceTwo.charAt(i - 1)))) {
                if(sequenceOne.charAt(j - 1) == sequenceTwo.charAt(i - 1)) {
                    comparator.insert(0, "|");
                    matchCount++;
                } else {
                    comparator.insert(0, "x");
                    misMatchCount++;
                }
                modifiedSequenceOne.insert(0, sequenceOne.charAt((j--) - 1));
                modifiedSequenceTwo.insert(0, sequenceTwo.charAt((i--) - 1));

                //  Else if gap in sequence two
            } else if(j > 0 && scoreMatrix[i][j] == scoreMatrix[i][j - 1] + ( i == 0 || i == lengthTwo ? initialGapPenalty : regularGapPenalty)) {
                comparator.insert(0, " ");
                modifiedSequenceOne.insert(0, sequenceOne.charAt((j--) - 1));
                modifiedSequenceTwo.insert(0, "_");
                gapCount++;

                //  Else if gap in sequence one
            } else if(i > 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j] + ( j == 0 || j == lengthOne ? initialGapPenalty : regularGapPenalty)) {
                comparator.insert(0, " ");
                modifiedSequenceOne.insert(0, "_");
                modifiedSequenceTwo.insert(0, sequenceTwo.charAt((i--) - 1));
                gapCount++;

                //  No previous cells could have lead here. Aborting.
            } else {
                throw new Exception("Error backtracking most optimal path. Aborting.");
            }
        }


        if(output.createNewFile()) {
//            System.out.println("Creating file '" + args[2] + "'.");
        } else {
//            System.out.println("Rewriting file '" + args[2] + "'.");
        }
        FileWriter fileWriter = new FileWriter(output);

        fileWriter.write(scoreMatrix[lengthTwo][lengthOne] + "\n");
        fileWriter.write(labelOne + "\n");
        fileWriter.write(modifiedSequenceOne + "\n");
        fileWriter.write(comparator + "\n");
        fileWriter.write(modifiedSequenceTwo + "\n");
        fileWriter.write(labelTwo.toString());

        fileWriter.close();
//        System.out.println("Saved results to " + output + ".\n");
//
//        System.out.println("Matches: " + matchCount);
//        System.out.println("Mismatches: " + misMatchCount);
//        System.out.println("Gaps: " + gapCount);
//        System.out.println("Gaps & Mismatches: " + (gapCount + misMatchCount));
//
//        System.out.println("\nScore: " + scoreMatrix[lengthTwo][lengthOne]);
    }

    static void parseDnaFromFile(File file, StringBuilder label, StringBuilder sequence) throws Exception {
        String line;
        boolean beginning = true;

        Scanner scanner = new Scanner(file);
        while(scanner.hasNextLine()) {
            line = scanner.nextLine().trim();
            if(line.charAt(0) == '>') {
                if(beginning) {
                    beginning = false;
                }
                label.append(line);
            } else {
                sequence.append(line);
                return;
            }
        }
        throw new Exception("Error while reading file: " + file + ".");
    }

    static int score(char base1, char base2) {
        if(base1 == base2) {
            return matchScore;
        } else {
            return mismatchPenalty;
        }
    }

    static StringBuilder trimToStartAndStop(StringBuilder stringBuilder) {
        StringBuilder newStringBuilder = new StringBuilder(stringBuilder);

        String currentCodon;
        int i = 0;
        for(; i < stringBuilder.length() - 3; i++) {
            currentCodon = stringBuilder.substring(i, i + 3);
            if(currentCodon.equals("ATG")) {
                break;
            }
        }
        if (i != 0) {
            newStringBuilder.delete(0, i);
        }
        int startIndex = i;
        for(; i < stringBuilder.length() - 3; i += 3) {
            currentCodon = stringBuilder.substring(i, i + 3);
            if(map.get(currentCodon).equals("O")) {
                break;
            }
        }
        newStringBuilder.delete(i - startIndex + 3, stringBuilder.length() - startIndex);
        return newStringBuilder;
    }

    static StringBuilder dnaToAa(StringBuilder dna) {
        StringBuilder dnaTrimmed = trimToStartAndStop(dna);
        StringBuilder aa = new StringBuilder();

        for(int i = 0; i < dnaTrimmed.length() - 2; i += 3) {
            aa.append(map.get(dnaTrimmed.substring(i, i + 3)));
        }

        return aa;
    }

    static void initHashmap() {
        map.put("TTT", "F");
        map.put("TTC", "F");
        map.put("TTA", "L");
        map.put("TTG", "L");
        map.put("CTT", "L");
        map.put("CTC", "L");
        map.put("CTA", "L");
        map.put("CTG", "L");
        map.put("ATT", "I");
        map.put("ATC", "I");
        map.put("ATA", "I");
        map.put("ATG", "M");
        map.put("GTT", "V");
        map.put("GTC", "V");
        map.put("GTA", "V");
        map.put("GTG", "V");
        map.put("TCT", "S");
        map.put("TCC", "S");
        map.put("TCA", "S");
        map.put("TCG", "S");
        map.put("CCT", "P");
        map.put("CCC", "P");
        map.put("CCA", "P");
        map.put("CCG", "P");
        map.put("ACT", "T");
        map.put("ACC", "T");
        map.put("ACA", "T");
        map.put("ACG", "T");
        map.put("GCT", "A");
        map.put("GCC", "A");
        map.put("GCA", "A");
        map.put("GCG", "A");
        map.put("TAT", "Y");
        map.put("TAC", "Y");
        map.put("TAA", "O");
        map.put("TAG", "O");
        map.put("TGA", "O");
        map.put("CAT", "H");
        map.put("CAC", "H");
        map.put("CAA", "Q");
        map.put("CAG", "Q");
        map.put("AAT", "N");
        map.put("AAC", "N");
        map.put("AAA", "K");
        map.put("AAG", "K");
        map.put("GAT", "D");
        map.put("GAC", "D");
        map.put("GAA", "E");
        map.put("GAG", "E");
        map.put("TGT", "C");
        map.put("TGC", "C");
        map.put("TGG", "W");
        map.put("CGT", "R");
        map.put("CGC", "R");
        map.put("CGA", "R");
        map.put("CGG", "R");
        map.put("AGT", "S");
        map.put("AGC", "S");
        map.put("AGA", "R");
        map.put("AGG", "R");
        map.put("GGT", "G");
        map.put("GGC", "G");
        map.put("GGA", "G");
        map.put("GGG", "G");
    }

}
