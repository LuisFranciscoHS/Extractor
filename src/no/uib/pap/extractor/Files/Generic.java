package no.uib.pap.extractor.Files;

import com.google.common.io.Files;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;

public class Generic {

    /**
     * Converts one file to the correct file format for PathwayMatcher
     *
     * @param args [0] input path+file [1] output path
     */
    public static void main(String args[]) {

//        splitFlexibleAndOne(args);

        convertToSimpleFormat(args);
    }

    private static void convertToSimpleFormat(String[] args) {
        FileWriter outputFlexible = null;
        if (args.length <= 0) {
            System.out.println("Need to specify the input file as argument.");
            System.exit(1);
        }
        try {
            List<String> lines = Files.readLines(new File(args[0]), Charset.defaultCharset());
            outputFlexible = new FileWriter(args[1] + "toMatchFlexible.tsv");

            int R = 0;
            for (String line : lines) {
                if (R == 0) {
                    R++;
                    continue;
                }
                if(line.contains("or")){
                    String[] parts = line.split("\t");
                    for(int I = 1; I < parts.length; I++){
                        outputFlexible.write(line.replace("\t", ";00000:").replace(" or ", ",00000:") + "\n");
                    }
                } else{
                    outputFlexible.write(line.replace("\t", ";00000:").replace(" and ", ",00000:") + "\n");
                }
            }
            outputFlexible.close();
            //outputOne.close();
        } catch (IOException e) {
            System.out.println("Could not read file: " + args[0]);
            System.exit(2);
        }
    }

    public static void splitFlexibleAndOne(String args[]) {

        FileWriter outputFlexible = null;
        FileWriter outputOne = null;

        if (args.length <= 0) {
            System.out.println("Need to specify the input file as argument.");
            System.exit(1);
        }
        try {
            List<String> lines = Files.readLines(new File(args[0]), Charset.defaultCharset());
            outputFlexible = new FileWriter(args[1] + "toMatchFlexible.tsv");
            outputOne = new FileWriter(args[1] + "toMatchOne.tsv");

            int R = 0;
            for (String line : lines) {
                if (R == 0) {
                    R++;
                    continue;
                }
                if(line.contains("or")){
                    outputOne.write(line.replace("\t", ";00000:").replace(" or ", ",00000:") + "\n");
                } else{
                    outputFlexible.write(line.replace("\t", ";00000:").replace(" and ", ",00000:") + "\n");
                }


            }
            outputFlexible.close();
            outputOne.close();
        } catch (IOException e) {
            System.out.println("Could not read file: " + args[0]);
            System.exit(2);
        }
    }

}
