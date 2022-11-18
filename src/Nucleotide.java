import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class Nucleotide {

    static String[] sampleNames;
    char[] nucleotides;  //reference, alternates
    int[][] numSamples; //sample index, nucleotide index
    double[][] frequencies;
    int[] depths;
    int position; //starts at index 1

    int codonNum;
    String referenceCodon = "";
    final static int MIN_DEPTH = 10;

    String gene;

    public Nucleotide(String data){
        String[] info = data.split("\t");
        position =Integer.parseInt(info[1]);
        initNucleotides(info);
        initDepthsAndFrequencies(info);
    }

    public Nucleotide(char[] nucleotides, int[][] numSamples, int[] depths, int position){
        this.nucleotides=nucleotides;
        this.numSamples = numSamples;
        this.depths=depths;
        this.position=position;
    }

    private void initDepthsAndFrequencies(String[] info){
        depths=new int[sampleNames.length];
        numSamples =new int[sampleNames.length][];
        frequencies = new double[sampleNames.length][];
        for (int i = 0; i < depths.length;i++) {
            String[] depthFrequency = info[i+9].split(":");
            depths[i]=Integer.parseInt(depthFrequency[1]);
            int[] num = new int[nucleotides.length];
            double[] frequency = new double[nucleotides.length];
            String[] stringFrequencies = depthFrequency[2].split(",");
            for (int j = 0; j < nucleotides.length;j++){
                num[j]=Integer.parseInt(stringFrequencies[j]);
                frequency[j]=(double)num[j]/depths[i];
            }
            numSamples[i]=num;
            frequencies[i]=frequency;
        }
    }

    private void initNucleotides(String[] info){
        char reference = info[3].charAt(0);
        String[] alternates = info[4].split(",");
        int numAlts = 0;
        for (int i = 0; i < alternates.length;i++){
            if (!alternates[i].equals("<*>")){
                numAlts++;
            }
        }
        nucleotides=new char[numAlts+1];
        nucleotides[0]=reference;
        for (int i = 0; i < alternates.length;i++){
            if (!alternates[i].equals("<*>")){
                nucleotides[i+1]=alternates[i].charAt(0);
            }
        }
    }

    public static void init(String header){

        initSampleNames(header);
    }

    private static void initSampleNames(String header){
        String[] names = header.split("\t");
        Nucleotide.sampleNames = new String[names.length-9];
        for (int i = 0; i < Nucleotide.sampleNames.length;i++){
            Nucleotide.sampleNames[i]=names[i+9];
        }
    }

    public double[] getReferenceFrequencies(){
        double[] referenceFrequencies = new double[sampleNames.length];
        for (int i = 0; i < referenceFrequencies.length;i++){
            if (depths[i]>=MIN_DEPTH) {
                referenceFrequencies[i] = frequencies[i][0];
            } else {
                referenceFrequencies[i]=-1;
            }
        }
        return referenceFrequencies;
    }

    public void addNucleotide(char nucleotide, int[] frequency){
        int index = indexOf(nucleotides,nucleotide);
        if (index==-1) {
            nucleotides = append(nucleotides, nucleotide);
            for (int i = 0; i < numSamples.length; i++) {
                numSamples[i] = append(numSamples[i], frequency[i]);
            }
        } else {
            for (int i = 0; i < numSamples.length; i++){
                numSamples[i][index]+=frequency[i];
            }
        }
    }

    public char[] getMutations(int sampleNum){
        if (depths[sampleNum]<MIN_DEPTH){
            return new char[]{nucleotides[0]};
        }
        int num = 0;
        for (int i = 0; i < frequencies[sampleNum].length;i++){
            if (frequencies[sampleNum][i]>=0.1){
                num++;
            }
        }
        char[] answer = new char[num];
        int index = 0;
        for (int i = 0; i < frequencies[sampleNum].length;i++){
            if (frequencies[sampleNum][i]>=0.1){
                answer[index]=nucleotides[i];
                index++;
            }
        }
        return answer;
    }

    private int indexOf(char[] array, char value){
        for (int i = 0; i < array.length;i++){
            if (array[i]==value){
                return i;
            }
        }
        return -1;
    }

    private int[] append(int[] array, int value){
        int[] answer = new int[array.length+1];
        for (int i = 0; i < array.length;i++){
            answer[i]=array[i];
        }
        answer[answer.length-1]=value;
        return answer;
    }

    private char[] append(char[] array, char value){
        char[] answer = new char[array.length+1];
        for (int i = 0; i < array.length;i++){
            answer[i]=array[i];
        }
        answer[answer.length-1]=value;
        return answer;
    }
}
