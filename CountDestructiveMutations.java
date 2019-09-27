/************************************************************           
 * MutPanning 									*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This scripts makes another quick run through	*
 * the maf file to count the number of protein destructive	*
 * mutaitons per gene. This mutaiton category will be used	*
 * as an indepdent test to detect tumor supressor geens. Pls*
 * note taht this is the only text in which insersion and	*
 * deletions are integrated into the mutational significance*
 * so that this test provides a non-redudant criterion. To	*
 * be robust against different annotations for nonsense		*
 * mutations in maf files, we determine this from the		*
 * reference sequence directly instead of from the maf file	*
 * Also associations of mutations with genes are determined	*
 * from the ref seq file instead of the maf file to ensure	*
 * consistency with the other tests and be indepednent of	*
 * different gene nomenclature systems.						*
 * 															*   
 *************************************************************/




import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

public class CountDestructiveMutations {

	static String[] entities=new String[0];
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};	
	static String[] chr2={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"};	

	static String file_align="";
	static String file_annotation="";
	static String file_maf="";
	static String file_entities="";
	static String file_genes="";
	static String file_out="";
	static String[] index_header_samples={"ID","Sample","Cohort"};
	static String[] index_header_maf={"Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode"};
	
	
	
	static ArrayList<Gene>[] genes_chr=new ArrayList[chr.length];
	static ArrayList<String>[] samples=new ArrayList[0];
	static ArrayList<Integer>[] samples_index=new ArrayList[0];
	static int[][][] n_indel=new int[chr.length][][];//[genes.size()][entities.length+1];
	static int[][][] n_indel_destructive=new int[chr.length][][];//[genes.size()][entities.length+1];
	static int[][][] mutations_snv_all=new int[chr.length][][];
	static int[][][] mutations_snv_nonsense=new int[chr.length][][];
	static int[][][] n_splice=new int[chr.length][][];//[genes.size()][entities.length+1];
	static int[][][] n_nonsense=new int[chr.length][][];//[genes.size()][entities.length+1];
	
	static String[] exclude={
			"PPIAL4A",//too long
			"PPIAL4B",
			"PPIAL4C",
			"FAM75A5",
			"FAM75A7",
			"ANKRD20A3",
			"FOXD4L2",
			"FOXD4L4",
			"FAM21B",
			"ASMT",//double gene
			"ASMTL",
			"CSF2RA",
			"DHRSX",
			"IL3RA",
			"IL3R",
			"IL9R",
			"SPRY3",
			"ZBED1"
		};
	static int[] chr_length={
			249250621,
			243199373,
			198022430,
			191154276,
			180915260,
			171115067,
			159138663,
			146364022,
			141213431,
			135534747,
			135006516,
			133851895,
			115169878,
			107349540,
			102531392,
			90354753,
			81195210,
			78077248,
			59128983,
			63025520,
			48129895,
			51304566,
			155270560,
			59373566
		};
	
	/*
	 * argument0:root file
	 * argument1: maf file
	 * argument2: sample file
	 * 
	 */
	
	public static void main(String[] args){
		
		file_annotation=args[3]+"AnnotationHg19/Annotation_chr";
		file_align=args[0]+"AlignHg19/AlignHg19Chr";
		file_maf=args[1];
		file_entities=args[2];
		file_genes=args[3]+"Exons_Hg19.txt";
		file_out=args[0]+"CountDestructive/";
		
		if(!new File(file_out).exists()){
			new File(file_out).mkdir();
		}
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
				try{
					FileInputStream in=new FileInputStream(file_entities);
					DataInputStream inn=new DataInputStream(in);
					BufferedReader input= new BufferedReader(new InputStreamReader(inn));
					int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
					String s="";
					ArrayList<String> aa =new ArrayList<String>();
					while((s=input.readLine())!=null){
						String e=s.split("	")[index_header[2]];
						if(!contains(e,aa)){
							aa.add(e);
						}
					}
					input.close();
					Collections.sort(aa);
					entities=new String[aa.size()];
					for (int i=0;i<aa.size();i++){
						entities[i]=aa.get(i);
					}
				}
				catch(Exception e){
					System.out.println(e);
				}		
				
				samples=new ArrayList[entities.length];
				samples_index=new ArrayList[entities.length];
		//System.out.println("X");
		try{
			
			//read samples and their index
			for (int i=0;i<samples.length;i++){
				samples[i]=new ArrayList<String>();
				samples_index[i]=new ArrayList<Integer>();
			}
			String s="";
			FileInputStream in=new FileInputStream(file_entities);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int ii=index(t[index_header[2]],entities);
				samples[ii].add(t[index_header[1]]);
				samples_index[ii].add(Integer.parseInt(t[index_header[0]]));
			}
			input.close();
			//System.out.println("XX");
			//read genes together with their coordinates
			for (int i=0;i<genes_chr.length;i++){
				genes_chr[i]=new ArrayList<Gene>();
			}
			
			in=new FileInputStream(file_genes);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				if(index(t[0],exclude)!=-1){
					continue;
				}
				//System.out.println(t[1]);
				int ii=index_gene(t[0],genes_chr[Integer.parseInt(t[1])-1]);
				if(ii==-1){
					genes_chr[Integer.parseInt(t[1])-1].add(new Gene(t[0],Integer.parseInt(t[2]),Integer.parseInt(t[3])));
				}
				else{
					genes_chr[Integer.parseInt(t[1])-1].get(ii).coord.add(new int[]{Integer.parseInt(t[2]),Integer.parseInt(t[3])});
				}
			}
			input.close();
			for (int i=0;i<genes_chr.length;i++){
				for (int j=0;j<genes_chr[i].size();j++){
					genes_chr[i].get(j).start=min(genes_chr[i].get(j).coord);
					genes_chr[i].get(j).end=max(genes_chr[i].get(j).coord);
				}
			}
			//System.out.println("XXX");
			
			//initalize counter arrays
			for (int i=0;i<chr.length;i++){
				n_indel[i]=new int[genes_chr[i].size()][entities.length+1];
				n_indel_destructive[i]=new int[genes_chr[i].size()][entities.length+1];
				mutations_snv_all[i]=new int[genes_chr[i].size()][entities.length+1];
				mutations_snv_nonsense[i]=new int[genes_chr[i].size()][entities.length+1];
				n_splice[i]=new int[genes_chr[i].size()][entities.length+1];
				n_nonsense[i]=new int[genes_chr[i].size()][entities.length+1];
				
			}
			System.out.println("XXX");
			
			Hashtable<String,Integer> table_samples=new Hashtable<String,Integer>();
			for (int i=0;i<samples.length;i++){
				for (int j=0;j<samples[i].size();j++){
					table_samples.put(samples[i].get(j), i);
				}
			}
			
			
			ArrayList<Integer> i_gene[][]=new ArrayList[chr.length][];
			for (int i=0;i<i_gene.length;i++){
				i_gene[i]=new ArrayList[1+chr_length[i]/10000];
				for (int j=0;j<i_gene[i].length;j++){
					i_gene[i][j]=new ArrayList<Integer>();
				}
				for (int j=0;j<genes_chr[i].size();j++){
					for (int k=genes_chr[i].get(j).start/10000;k<=genes_chr[i].get(j).end/10000;k++){
						i_gene[i][k].add(j);
					}
				}
				
			}
			
			
			//go through maf file and extract count information 
			in=new FileInputStream(file_maf);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header_m=index_header(input.readLine().split("	"),index_header_maf);
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				
				int jj=table_samples.get(t[index_header_m[10]]).intValue();//index(,samples);
				int chr_index=index(t[index_header_m[1]],chr,chr2);
				if(chr_index==-1){
					continue;
				}
				ArrayList<Integer> gene_index=index(Integer.parseInt(t[index_header_m[2]]),genes_chr[chr_index],i_gene[chr_index][Integer.parseInt(t[index_header_m[2]])/10000]);
				String ref=t[index_header_m[7]].toUpperCase();
				String tumor=t[index_header_m[9]].toUpperCase();
				String type=t[index_header_m[5]];
				
				if(tumor.equals(ref)||tumor.equals("")){
					tumor=t[index_header_m[8]].toUpperCase();
				}
				
				for (int i=0;i<gene_index.size();i++){
					if(ref.equals("-")){//Inserstion
						n_indel[chr_index][gene_index.get(i)][jj]++;
						n_indel[chr_index][gene_index.get(i)][entities.length]++;
						if(tumor.length()%3!=0){//frameshift
							n_indel_destructive[chr_index][gene_index.get(i)][jj]++;
							n_indel_destructive[chr_index][gene_index.get(i)][entities.length]++;
						}
					}
					else if(tumor.equals("-")){//Deletion
						n_indel[chr_index][gene_index.get(i)][jj]++;
						n_indel[chr_index][gene_index.get(i)][entities.length]++;
						if(ref.length()%3!=0){//frameshift
							n_indel_destructive[chr_index][gene_index.get(i)][jj]++;
							n_indel_destructive[chr_index][gene_index.get(i)][entities.length]++;
						}
					}
				}
				
				if(type.toLowerCase().contains("splice")){
					
					int index_gene=index_gene(t[index_header_m[0]],genes_chr[chr_index]);
					if(index_gene!=-1){
						n_splice[chr_index][index_gene][jj]++;
						n_splice[chr_index][index_gene][entities.length]++;
					}
					
				}
				if(type.toLowerCase().contains("nonsense")){
					
					int index_gene=index_gene(t[index_header_m[0]],genes_chr[chr_index]);
					//System.out.println(type+"	"+index_gene);
					if(index_gene!=-1){
						n_nonsense[chr_index][index_gene][jj]++;
						n_nonsense[chr_index][index_gene][entities.length]++;
					}
					
				}
			}
			input.close();
			System.out.println("XXX");
			
			
			//go through the reference sequence of each chromosome and
			//look for nonsynonymous mutations based on the reference seq
			for (int i=0;i<chr.length;i++){
				System.out.println(chr[i]);
				run_subthread(i);
			}
			/*
			Subthread[] threads=new Subthread[chr.length];
			for (int i=0;i<threads.length;i++){
				threads[i]=new Subthread();
				threads[i].c=i;
				threads[i].start();
			}

			//wait for all subthreads
			boolean all_done=false;
			do{
				Thread.sleep(3000);
				all_done=true;
				for (int i=0;i<threads.length;i++){
					if(!threads[i].done){
						all_done=false;
						break;
					}
				}
			}while(!all_done);
			*/
			//output for each gene of the number of all mutations  and the number of destructive mutaitons
			for (int k=0;k<entities.length;k++){
				FileWriter out=new FileWriter(file_out+entities[k]+".txt");
				BufferedWriter output= new BufferedWriter(out);				
				for (int i=0;i<genes_chr.length;i++){
					for (int j=0;j<genes_chr[i].size();j++){
						//output.write(genes_chr[i].get(j).name+"	"+(n_indel[i][j][k]+mutations_snv_all[i][j][k])+"	"+(n_indel_destructive[i][j][k]+mutations_snv_nonsense[i][j][k]));
						if(mutations_snv_nonsense[i][j][k]<n_nonsense[i][j][k]){
							mutations_snv_nonsense[i][j][k]=n_nonsense[i][j][k];
						}
						if(mutations_snv_nonsense[i][j][k]>mutations_snv_all[i][j][k]){
							mutations_snv_all[i][j][k]=mutations_snv_nonsense[i][j][k];
						}
								
						
						output.write(genes_chr[i].get(j).name+"	"+(n_indel[i][j][k]+mutations_snv_all[i][j][k]+n_splice[i][j][k])+"	"+(n_indel[i][j][k]+mutations_snv_nonsense[i][j][k]+n_splice[i][j][k]));
						output.newLine();
					}
				}
				output.close();
			}
			FileWriter out=new FileWriter(file_out+"PanCancer"+".txt");
			BufferedWriter output= new BufferedWriter(out);				
			for (int i=0;i<genes_chr.length;i++){
				for (int j=0;j<genes_chr[i].size();j++){
					//output.write(genes_chr[i].get(j).name+"	"+(n_indel[i][j][entities.length]+mutations_snv_all[i][j][entities.length])+"	"+(n_indel_destructive[i][j][entities.length]+mutations_snv_nonsense[i][j][entities.length]));
					if(mutations_snv_nonsense[i][j][entities.length]<n_nonsense[i][j][entities.length]){
						mutations_snv_nonsense[i][j][entities.length]=n_nonsense[i][j][entities.length];
					}
					if(mutations_snv_nonsense[i][j][entities.length]>mutations_snv_all[i][j][entities.length]){
						mutations_snv_all[i][j][entities.length]=mutations_snv_nonsense[i][j][entities.length];
					}
					output.write(genes_chr[i].get(j).name+"	"+(n_indel[i][j][entities.length]+mutations_snv_all[i][j][entities.length]+n_splice[i][j][entities.length])+"	"+(n_indel[i][j][entities.length]+mutations_snv_nonsense[i][j][entities.length]+n_splice[i][j][entities.length]));
					
					output.newLine();
				}
			}
			output.close();
			
			
		}
		catch(Exception e){
			
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
		}
	}
	
	public static int[] index_header(String[] header, String[] ideal_header){
		int[] indices=new int[ideal_header.length];
		for (int i=0;i<ideal_header.length;i++){
			int index=-1;
			for (int j=0;j<header.length;j++){
				if(header[j].equals(ideal_header[i])){
					index=j;
					break;
				}
			}
			indices[i]=index;
		}
		return indices;
	}
	
	public static ArrayList<Integer> index(int pos,ArrayList<Gene> genes,ArrayList<Integer> a){
		ArrayList<Integer> index=new ArrayList<Integer>();
		for (int ii=0;ii<a.size();ii++){
			if(genes.get(a.get(ii)).contains(pos)){
				index.add(a.get(ii));
			}
		}
		return index;
	}
	
	public static boolean contains (String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	public static int index (String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	public static int index (String s, ArrayList<String>[] t){
		for (int i=0;i<t.length;i++){
			for (int j=0;j<t[i].size();j++){
				if(t[i].get(j).equals(s)){
					return i;
				}
			}
		}
		
		return -1;
	}
	
	public static boolean contains (String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return true;
			}
		}
		return false;
	}
	
	public static int index (String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	public static int index (String s, String[] t1, String[] t2){
		int i1=index(s,t1);
		if(i1==-1){
			return index(s,t2);
		}
		return i1;
	}
	
	public static int min(ArrayList<int[]> a){
		int min=1000000000;
		for (int i=0;i<a.size();i++){
			if(a.get(i)[0]<min){
				min=a.get(i)[0];
			}
		}
		return min;
	}
	
	public static int max(ArrayList<int[]> a){
		int max=-1000000000;
		for (int i=0;i<a.size();i++){
			if(a.get(i)[1]>max){
				max=a.get(i)[1];
			}
		}
		return max;
	}
	
	
	//go through the reference sequence of each gene and 
	//search for nonsynonymous mutations
	
	
	
	/*
	public static void run_subthread(int c){
		System.out.println("start "+chr[c]);
		try{
			Hashtable<Integer,Integer> table_samples_index=new Hashtable<Integer,Integer>();
			for (int i=0;i<samples_index.length;i++){
				for (int j=0;j<samples_index[i].size();j++){
					table_samples_index.put(samples_index[i].get(j),i);
				}
			}
			
			ArrayList<Integer> position=new ArrayList<Integer>();
	//		ArrayList<String> nucl=new ArrayList<String>();
	//		ArrayList<Double> coverage=new ArrayList<Double>();
			ArrayList<int[][]> count=new ArrayList<int[][]>();
			
			FileInputStream in2=new FileInputStream(file_align+chr[c]+".txt");
			DataInputStream inn2=new DataInputStream(in2);
			BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
			
			for (int nn=0;nn<genes_chr[c].size();nn++){
				System.out.println(nn+"/"+genes_chr[c].size());
				int ii=0;
				if(position.size()>0){
					while(ii<position.size()&&position.get(ii)<genes_chr[c].get(nn).start){
						ii++;
					}
				}
				if(ii>0){
					for (int i=ii-1;i>=0;i--){
						//nucl.remove(i);
						position.remove(i);
						//coverage.remove(i);
						count.remove(i);
						
					}
				}
				
				String s2="";
				while ((s2=input2.readLine())!=null){
					
					
					String[] t2=s2.split("	");
					
					
					String ref_as="";
					String mut_as1="";
					String mut_as2="";
					String mut_as3="";
					if(t2.length>6){
						ref_as=t2[6];
						mut_as1=t2[8];
						mut_as2=t2[9];
						mut_as3=t2[10];
					}
					int[] nonsense=new int[3];
					if(ref_as.equals("*")||mut_as1.equals("*")){
						nonsense[0]=1;
					}
					if(ref_as.equals("*")||mut_as2.equals("*")){
						nonsense[1]=1;
					}
					if(ref_as.equals("*")||mut_as3.equals("*")){
						nonsense[2]=1;
					}
					
					
					int[] index1=new int[0];
					if(t2.length>3&&!t2[3].equals("")){
						index1=integer(t2[3].split(";"));
					}
							
					int[] index2=new int[0];
					if(t2.length>4&&!t2[4].equals("")){
						index2=integer(t2[4].split(";"));
					}
					int[] index3=new int[0];
					if(t2.length>5&&!t2[5].equals("")){
						index3=integer(t2[5].split(";"));
					}
					int[][] indices=new int[][]{index1,index2,index3};
					
					
					int[][] count_local=new int[entities.length+1][2];
					for (int k=0;k<indices.length;k++){
						
						for (int l=0;l<indices[k].length;l++){
							count_local[table_samples_index.get(indices[k][l])][nonsense[k]]++;
						}
						//for (int i=0;i<samples.length;i++){
							//count_local[i][nonsense[k]]+=overlap(indices[k],samples_index[i]);
						//}
						count_local[count_local.length-1][nonsense[k]]=indices[k].length;
					}
					count.add(count_local);
					
					
					position.add(Integer.parseInt(t2[0]));
					//nucl.add(t2[1]);
					//coverage.add(Double.parseDouble(t2[2]));
					
					
					if(Integer.parseInt(t2[0])>=genes_chr[c].get(nn).end){
						break;
					}
				}
				
				
				for (int i=0;i<position.size();i++){
					
					if(genes_chr[c].get(nn).contains(position.get(i))){
						for (int k=0;k<mutations_snv_all[c][nn].length;k++){
							mutations_snv_all[c][nn][k]+=count.get(i)[k][0]+count.get(i)[k][1];
							mutations_snv_nonsense[c][nn][k]+=count.get(i)[k][1];
						}
					}
				}
				
				
			}
			input2.close();
			
		}
		catch(Exception e){
			
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
		}
		System.out.println("done "+chr[c]);

	}
	*/
	
	public static void run_subthread(int c){
		System.out.println("start "+chr[c]);
		try{
			Hashtable<Integer,Integer> table_samples_index=new Hashtable<Integer,Integer>();
			for (int i=0;i<samples_index.length;i++){
				for (int j=0;j<samples_index[i].size();j++){
					table_samples_index.put(samples_index[i].get(j),i);
				}
			}
			
			ArrayList<Integer> position=new ArrayList<Integer>();
			ArrayList<Integer> nonsense=new ArrayList<Integer>();
			ArrayList<Integer> entity=new ArrayList<Integer>();
			
	//		ArrayList<String> nucl=new ArrayList<String>();
	//		ArrayList<Double> coverage=new ArrayList<Double>();
	//		ArrayList<int[][]> count=new ArrayList<int[][]>();
			
			FileInputStream in2=new FileInputStream(file_align+chr[c]+".txt");
			DataInputStream inn2=new DataInputStream(in2);
			BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
			
			FileInputStream in3=new FileInputStream(file_annotation+chr[c]+".txt");
			DataInputStream inn3=new DataInputStream(in3);
			BufferedReader input3= new BufferedReader(new InputStreamReader(inn3));
			
			
			for (int nn=0;nn<genes_chr[c].size();nn++){
				//System.out.println(nn+"/"+genes_chr[c].size());
				int ii=0;
				if(position.size()>0){
					while(ii<position.size()&&position.get(ii)<genes_chr[c].get(nn).start){
						ii++;
					}
				}
				if(ii>0){
					for (int i=ii-1;i>=0;i--){
						//nucl.remove(i);
						position.remove(i);
						nonsense.remove(i);
						entity.remove(i);
						//coverage.remove(i);
						//count.remove(i);
						
					}
				}
				
				String s2="";
				while ((s2=input2.readLine())!=null){
					
					
					String[] t2=s2.split("	");
					String[] t3=input3.readLine().split("	");
					
					
					String ref_as="";
					String mut_as1="";
					String mut_as2="";
					String mut_as3="";
					if(t3.length>3){
						ref_as=t3[3];
						mut_as1=t3[5];
						mut_as2=t3[6];
						mut_as3=t3[7];
					}
					
					
					/*
					int[] nonsense=new int[3];
					if(ref_as.equals("*")||mut_as1.equals("*")){
						nonsense[0]=1;
					}
					if(ref_as.equals("*")||mut_as2.equals("*")){
						nonsense[1]=1;
					}
					if(ref_as.equals("*")||mut_as3.equals("*")){
						nonsense[2]=1;
					}*/
					
					int nonsense1=0;
					int nonsense2=0;
					int nonsense3=0;
					if(ref_as.equals("*")||mut_as1.equals("*")){
						nonsense1=1;
					}
					if(ref_as.equals("*")||mut_as2.equals("*")){
						nonsense2=1;
					}
					if(ref_as.equals("*")||mut_as3.equals("*")){
						nonsense3=1;
					}
					
					if(t2.length>0&&!t2[0].equals("")){
						String[] tt=t2[0].split(";");
						for (int i=0;i<tt.length;i++){
							entity.add(table_samples_index.get(Integer.parseInt(tt[i])));
							nonsense.add(nonsense1);
							position.add(Integer.parseInt(t3[0]));
						}
					}
					if(t2.length>1&&!t2[1].equals("")){
						String[] tt=t2[1].split(";");
						for (int i=0;i<tt.length;i++){
							entity.add(table_samples_index.get(Integer.parseInt(tt[i])));
							nonsense.add(nonsense2);
							position.add(Integer.parseInt(t3[0]));
						}
					}
					if(t2.length>2&&!t2[2].equals("")){
						String[] tt=t2[2].split(";");
						for (int i=0;i<tt.length;i++){
							entity.add(table_samples_index.get(Integer.parseInt(tt[i])));
							nonsense.add(nonsense3);
							position.add(Integer.parseInt(t3[0]));
						}
					}
					
					/*
					String[] tt1=t2[3].split(";");
					String[] tt2=t2[4].split(";");
					String[] tt3=t2[5].split(";");
					
					
					
					
					int[] index1=new int[0];
					if(t2.length>3&&!t2[3].equals("")){
						index1=integer(t2[3].split(";"));
					}
							
					int[] index2=new int[0];
					if(t2.length>4&&!t2[4].equals("")){
						index2=integer(t2[4].split(";"));
					}
					int[] index3=new int[0];
					if(t2.length>5&&!t2[5].equals("")){
						index3=integer(t2[5].split(";"));
					}
					int[][] indices=new int[][]{index1,index2,index3};
					
					
					int[][] count_local=new int[entities.length+1][2];
					for (int k=0;k<indices.length;k++){
						
						for (int l=0;l<indices[k].length;l++){
							count_local[table_samples_index.get(indices[k][l])][nonsense[k]]++;
						}
						//for (int i=0;i<samples.length;i++){
							//count_local[i][nonsense[k]]+=overlap(indices[k],samples_index[i]);
						//}
						count_local[count_local.length-1][nonsense[k]]=indices[k].length;
					}
					count.add(count_local);
					
					
					position.add(Integer.parseInt(t2[0]));
					//nucl.add(t2[1]);
					//coverage.add(Double.parseDouble(t2[2]));
					*/
					
					if(Integer.parseInt(t3[0])>=genes_chr[c].get(nn).end){
						break;
					}
				}
				
				
				for (int i=0;i<position.size();i++){
					if(genes_chr[c].get(nn).contains(position.get(i))){
						mutations_snv_all[c][nn][entity.get(i)]++;
						if(nonsense.get(i)==1){
							mutations_snv_nonsense[c][nn][entity.get(i)]++;
						}
						mutations_snv_all[c][nn][entities.length]++;
						if(nonsense.get(i)==1){
							mutations_snv_nonsense[c][nn][entities.length]++;
						}
						
						/*
						for (int k=0;k<mutations_snv_all[c][nn].length;k++){
							mutations_snv_all[c][nn][k]+=count.get(i)[k][0]+count.get(i)[k][1];
							mutations_snv_nonsense[c][nn][k]+=count.get(i)[k][1];
						}*/
					}
				}
				
				
			}
			input2.close();
			input3.close();
		}
		catch(Exception e){
			
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
		}
		System.out.println("done "+chr[c]);

	}

	public static int[] integer(String[] s){
		int[] a=new int[s.length];
		for (int i=0;i<a.length;i++){
			a[i]=Integer.parseInt(s[i]);
		}
		return a;
	}
	private static class Gene{
		String name="";
		ArrayList<int[]> coord=new ArrayList<int[]>();
		int start=-1;
		int end=-1;
		public Gene(String name, int start, int end){
			this.name=name;
			coord.add(new int[]{start,end});
		}
		public boolean contains (int pos){
			for (int i=0;i<coord.size();i++){
				if(coord.get(i)[0]<=pos&&pos<=coord.get(i)[1]){
					return true;
				}
			}
			return false;
		}
	}
	public static int index_gene(String s, ArrayList<Gene> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).name.equals(s)){
				return i;
			}
		}
		return -1;
	}
	public static int overlap(int[] a, ArrayList<Integer> b){
		int n=0;
		for (int i=0;i<a.length;i++){
			for (int j=0;j<b.size();j++){
				if(a[i]==b.get(j)){
					n++;
				}
			}
		}
		return n;
	}
}
