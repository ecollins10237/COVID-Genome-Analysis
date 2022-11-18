import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Scanner;

public class Main {

    public static void main(String[] args) throws IOException {
        String[] genes = readAnnotation();
        String filePath="COVID_COMPLETECURES_d3000_tDPAD_v2.vcf[79]";
        Nucleotide[] genome = readFile(filePath);
        generateSynonymousMatrix(genome,genes,filePath);
        generateRAFMatrix(genome,filePath);
    }

    public static Nucleotide[]  readFile(String filePath) throws IOException {
        ArrayList<Nucleotide> nucleotides = new ArrayList<>();
        File f = new File(filePath);
        Scanner scanner = new Scanner(f);
        String line = scanner.nextLine();
        while (!line.contains("#CHROM")){
            line=scanner.nextLine();
        }
        Codon.init(line);
        Nucleotide.init(line);
        ArrayList<String> INDELS = new ArrayList<>();
        while (scanner.hasNextLine()){
            line=scanner.nextLine();
            if (line.split("\t")[3].length()==1) {
                nucleotides.add(new Nucleotide(line));
            } else {
                INDELS.add(line);
            }
        }
        Nucleotide[] genome = new Nucleotide[nucleotides.size()];
        for (int i = 0; i < genome.length;i++){
            genome[i]=nucleotides.get(i);
        }
        return genome;
    }

    public static void generateSynonymousMatrix(Nucleotide[] nucleotides, String[] genes, String fileName) throws IOException {
        Codon[] codons = new Codon[nucleotides.length/3];
        for (int i = 0; i < nucleotides.length-3;i+=3){
            codons[i/3]=new Codon(nucleotides[i],nucleotides[i+1],nucleotides[i+2]);
        }
        fileName=fileName.substring(0,fileName.indexOf(".vcf"))+"SN";
        FileWriter fileWriter = new FileWriter(fileName+".csv");
        String header = " ,"+Arrays.toString(Nucleotide.sampleNames).substring(1,Arrays.toString(Nucleotide.sampleNames).length()-1);
        fileWriter.write(header+"\n");
        for (int i = 0; i < codons.length;i++){
            String line = genes[codons[i].n1.position-1]+"_"+Codon.codonTable.get(codons[i].reference)+"_"+(i+1)+",";
            for (int j = 0; j < Codon.sampleNames.length;j++){
                String aminoAcidList = "";
                ArrayList<String> aminoAcids = codons[i].getMutatedAminoAcids(j);
                for (int k = 0; k < aminoAcids.size();k++){
                    aminoAcidList+=aminoAcids.get(k)+":";
                }
                line+=aminoAcidList.substring(0,Math.max(aminoAcidList.length()-1,0))+",";
            }
            fileWriter.write(line.substring(0,line.length()-1)+"\n");
        }
        fileWriter.close();
    }

    public static void generateRAFMatrix(Nucleotide[] nucleotides, String fileName) throws IOException {
        fileName=fileName.substring(0,fileName.indexOf(".vcf"))+"RAF";
        FileWriter fileWriter = new FileWriter(fileName+".csv");
        String header = "POS,";
        for (int i = 0; i < Nucleotide.sampleNames.length;i++){
            header+=Nucleotide.sampleNames[i]+",";
        }
        fileWriter.write(header.substring(0,header.length()-1)+"\n");
        for (int i = 0; i < nucleotides.length;i++){
            String line = nucleotides[i].position +",";
            double[] frequencies = nucleotides[i].getReferenceFrequencies(); //depth >= 10
            for (int j = 0; j < frequencies.length;j++){
                line+=frequencies[j]==-1? "NaN,":(frequencies[j]+",");
            }
            fileWriter.write(line.substring(0,line.length()-1)+"\n");
        }
        fileWriter.close();
    }

    private static String[] readAnnotation() throws FileNotFoundException {
        ArrayList<String> genes = new ArrayList<>();
        Scanner file = new Scanner(new File("COVID_Annotation.csv"));
        file.nextLine();
        while (file.hasNextLine()) {
            String line = file.nextLine();
            if (line.contains(",")) {
                String[] info = line.split(",");
                genes.add(info[3]);
            }
        }
        return genes.toArray(new String[0]);
    }
}