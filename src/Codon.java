import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

public class Codon {

    static String[] sampleNames;

    Nucleotide n1,n2,n3;
    String reference;
    String[][] mutations;
    static Hashtable<String,String> codonTable = new Hashtable<>();

    public Codon(Nucleotide n1, Nucleotide n2, Nucleotide n3){
        this.n1=n1;
        this.n2=n2;
        this.n3=n3;
        reference=""+n1.nucleotides[0]+n2.nucleotides[0]+n3.nucleotides[0];
        mutations=new String[sampleNames.length][];
        for (int i = 0; i < sampleNames.length;i++){
            char[] n1M = n1.getMutations(i);
            char[] n2M = n2.getMutations(i);
            char[] n3M = n3.getMutations(i);
            String[] m = new String[n1M.length*n2M.length*n3M.length];
            int counter = 0;
            for (int a = 0; a<n1M.length;a++){
                for (int b = 0; b < n2M.length;b++){
                    for (int c = 0; c < n3M.length; c++){
                        m[counter] = ""+n1M[a]+n2M[b]+n3M[c];
                        counter++;
                    }
                }
            }
            if (m.length==0){
                mutations[i]=new String[]{reference};
            } else {
                mutations[i] = m;
            }
        }
    }

    public String getSynonymous(int sampleNum){
        if (mutations[sampleNum].length==1&&mutations[sampleNum][0].equals(reference)){
            return "";
        }
        String referenceAmino = codonTable.get(reference);
        for (int i = 0; i < mutations[sampleNum].length;i++){
            if (!codonTable.get(mutations[sampleNum][i]).equals(referenceAmino)){
                return "N";
            }
        }
        return "S";
    }

    public ArrayList<String> getMutatedAminoAcids(int sampleNum){
        ArrayList<String> answer = new ArrayList<>();
        for (int i = 0; i < mutations[sampleNum].length;i++){
            if (!answer.contains(codonTable.get(mutations[sampleNum][i]))){
                answer.add(codonTable.get(mutations[sampleNum][i]));
            }
        }
        return answer;
    }

    public static void init(String header){
        initCodonTable();
        initSampleNames(header);
    }

    private static void initSampleNames(String header){
        String[] names = header.split("\t");
        Codon.sampleNames = new String[names.length-9];
        for (int i = 0; i < Codon.sampleNames.length;i++){
            Codon.sampleNames[i]=names[i+9];
        }
    }

    private static void initCodonTable(){
        codonTable=new Hashtable<>();
        final String[] aminoNames = {"F","L","S","Y","*","C","*","W","L","P","H","Q","R","I","M","T","N","K","S","R","V","A","D","E","G"};
        final int[] codonsPerAmino = {2,2,4,2,2,2,1,1,4,4,2,2,4,3,1,4,2,2,2,2,4,4,2,2,4};
        final char[] bases = {'T','C','A','G'};
        int count = 0;
        int index = 0;
        for (int i = 0; i < 4;i++){
            for (int j = 0; j < 4;j++){
                for (int k = 0; k < 4;k++){
                    codonTable.put(""+bases[i]+bases[j]+bases[k],aminoNames[index]);
                    count++;
                    if (count>=codonsPerAmino[index]){
                        index++;
                        count=0;
                    }
                }
            }
        }
    }
}
