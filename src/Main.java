import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class Main {

    public static void main(String[] args) throws IOException {
        String filePath="COVID_COMPLETECURES_d3000_tDPAD_v2.vcf[79]";
        Nucleotide[] genome = readData(filePath);
        generateRAFMatrix(genome,filePath);
        genome = readSNData(filePath);
        generateSNSynonymousMatrix(genome,filePath);
        generateAminoAcidSynonymousMatrix(genome,filePath);
    }

    public static Nucleotide[] readData(String filePath) throws IOException {
        ArrayList<Nucleotide> nucleotides = new ArrayList<>();
        Scanner data = new Scanner(new File(filePath));
        String line = data.nextLine();
        while (!line.contains("#CHROM")){
            line=data.nextLine();
        }
        Codon.init(line);
        Nucleotide.init(line);
        while (data.hasNextLine()){
            line=data.nextLine();
            if (line.split("\t")[3].length()==1) {
                nucleotides.add(new Nucleotide(line));
            }
        }
        Nucleotide[] genome = new Nucleotide[nucleotides.size()];
        for (int i = 0; i < genome.length;i++){
            genome[i]=nucleotides.get(i);
        }
        return genome;
    }

    public static Nucleotide[] readSNData(String filePath) throws IOException {
        ArrayList<Nucleotide> nucleotides = new ArrayList<>();
        Scanner data = new Scanner(new File(filePath));
        Scanner annotations = new Scanner(new File("COVID_Annotation.csv"));
        annotations.nextLine();
        annotations.nextLine();
        String line = data.nextLine();
        while (!line.contains("#CHROM")){
            line=data.nextLine();
        }
        Codon.init(line);
        Nucleotide.init(line);
        ArrayList<Nucleotide> toCopy = new ArrayList<>();
        while (data.hasNextLine()){
            line=data.nextLine();
            if (!line.contains("INDEL")) {
                Nucleotide nucleotide = new Nucleotide(line);
                String[] info = annotations.nextLine().split(",");
                if (info[3].contains("/")) {
                    nucleotide.gene = info[3].substring(0, info[3].indexOf("/"));
                    nucleotide.codonNum = Integer.parseInt(info[6].substring(0,info[6].indexOf("/")));
                    nucleotide.referenceCodon = info[7].substring(0,info[7].indexOf("/"));
                    nucleotides.add(nucleotide);
                    toCopy.add(new Nucleotide(line));
                    toCopy.get(toCopy.size() - 1).gene = info[3].substring(info[3].indexOf("/")+1);
                    toCopy.get(toCopy.size() - 1).codonNum = Integer.parseInt(info[6].substring(info[6].indexOf("/")+1));
                    toCopy.get(toCopy.size() - 1).referenceCodon = info[7].substring(info[7].indexOf("/")+1);
                } else {
                    for (int i = 0; i < toCopy.size();i++){
                        nucleotides.add(toCopy.get(i));
                    }
                    toCopy.clear();
                    nucleotide.gene = info[3];
                    if (info.length >= 6) {
                        nucleotide.codonNum = Integer.parseInt(info[6]);
                        nucleotide.referenceCodon = info[7];
                    }
                    nucleotides.add(nucleotide);
                }
                if (annotations.hasNextLine()){
                    annotations.nextLine();
                }
            }
        }
        Nucleotide[] genome = new Nucleotide[nucleotides.size()];
        for (int i = 0; i < genome.length;i++){
            genome[i]=nucleotides.get(i);
        }
        return genome;
    }

    public static void generateSNSynonymousMatrix(Nucleotide[] nucleotides, String fileName) throws IOException {
        System.out.print("Generating Synonomous Matrix S and N . . . ");
        long startTime =  System.currentTimeMillis();
        ArrayList<Codon> codons = new ArrayList<>();
        for (int i = 0; i < nucleotides.length-3;i+=3){
            if (nucleotides[i].referenceCodon.length()!=0) {
                codons.add(new Codon(nucleotides[i], nucleotides[i + 1], nucleotides[i + 2]));
            } else {
                i-=2;
            }
        }
        fileName=fileName.substring(0,fileName.indexOf(".vcf"))+"SN_"+Nucleotide.MIN_DEPTH+"_"+Nucleotide.MIN_PERCENT;;
        FileWriter fileWriter = new FileWriter(fileName+".csv");
        String header = " ,"+Arrays.toString(Nucleotide.sampleNames).substring(1,Arrays.toString(Nucleotide.sampleNames).length()-1);
        fileWriter.write(header+"\n");
        for (int i = 0; i < codons.size();i++){
            String line = codons.get(i).n1.gene+"_"+Codon.codonTable.get(codons.get(i).reference)+"_"+codons.get(i).n1.codonNum+",";
            for (int j = 0; j < Codon.sampleNames.length;j++){
                line+=codons.get(i).getSynonymous(j)+",";
            }
            fileWriter.write(line.substring(0,line.length()-1)+"\n");
        }
        fileWriter.close();
        System.out.println(System.currentTimeMillis()-startTime+" milliseconds");
    }

    public static void generateAminoAcidSynonymousMatrix(Nucleotide[] nucleotides, String fileName) throws IOException {
        System.out.print("Generating Synonomous Matrix with Amino Acid Names . . . ");
        long startTime = System.currentTimeMillis();
        ArrayList<Codon> codons = new ArrayList<>();
        for (int i = 0; i < nucleotides.length-3;i+=3){
            if (nucleotides[i].referenceCodon.length()!=0) {
                codons.add(new Codon(nucleotides[i], nucleotides[i + 1], nucleotides[i + 2]));
            } else {
                i-=2;
            }
        }
        fileName=fileName.substring(0,fileName.indexOf(".vcf"))+"SNAminoAcid_"+Nucleotide.MIN_DEPTH+"_"+Nucleotide.MIN_PERCENT;;
        FileWriter fileWriter = new FileWriter(fileName+".csv");
        String header = " ,"+Arrays.toString(Nucleotide.sampleNames).substring(1,Arrays.toString(Nucleotide.sampleNames).length()-1);
        fileWriter.write(header+"\n");
        for (int i = 0; i < codons.size();i++){
            String line = codons.get(i).n1.gene+"_"+Codon.codonTable.get(codons.get(i).reference)+"_"+codons.get(i).n1.codonNum+",";
            for (int j = 0; j < Codon.sampleNames.length;j++){
                String synonomous = codons.get(i).getSynonymous(j);
                if (!synonomous.equals("")){
                    String acids="";
                    ArrayList<String> mutatedAcids = codons.get(i).getMutatedAminoAcids(j);
                    for (int k = 0; k < mutatedAcids.size();k++){
                        acids+=mutatedAcids.get(k)+":";
                    }
                    line+=acids.substring(0,acids.length()-1);
                }
                line+=",";
            }
            fileWriter.write(line.substring(0,line.length()-1)+"\n");
        }
        fileWriter.close();
        System.out.println(System.currentTimeMillis()-startTime+" milliseconds");
    }

    public static void generateRAFMatrix(Nucleotide[] nucleotides, String fileName) throws IOException {
        System.out.print("Generating RAF Matrix . . . ");
        long startTime = System.currentTimeMillis();
        fileName=fileName.substring(0,fileName.indexOf(".vcf"))+"RAF_"+Nucleotide.MIN_DEPTH+"_"+Nucleotide.MIN_PERCENT;
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
        System.out.println(System.currentTimeMillis()-startTime+" milliseconds");
    }
}