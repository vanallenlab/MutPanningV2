/************************************************************           
 * MutPanning 												*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: Another quick auxiliary script which sums		*
 * the mutation counts of the previous steps which were		*
 * generated in each chr separately into a single file and	*
 * separately for syn and nonsynonymous counts, as needed	*
 * for the CBASE model										*
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


public class ReformatCBASE{
	
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};	
	static String[] entity_name=new String[0];
	
	static String[] index_header_samples={"ID","Sample","Cohort"};
	static String[] index_header_maf={"Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode"};
	
	
	static String file_count="";
	static String file_maf="";
	static String file_genes="";
	static String file_samples="";
	
	static String file_out1="";
	static String file_out2="";
	static String file_out3="";
	
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
			"ASMT",//exists twice 
			"ASMTL",
			"CSF2RA",
			"DHRSX",
			"IL3RA",
			"IL3R",
			"IL9R",
			"SPRY3",
			"ZBED1"
		};
	
	/*
	 * argument0: root file 
	 * argument1: maf file
	 * argument2: sample count file
	 */
	public static void main(String[] args){
		
		file_count=args[0]+"CBASE/CountsRaw/Count";
		file_maf=args[1];
		file_genes=args[3]+"Exons_Hg19.txt";
		file_samples=args[2];
		
		file_out1=args[0]+"CBASE/Counts/Count";
		file_out2=args[0]+"CBASE/Counts/CountSilent";
		file_out3=args[0]+"CBASE/CountsChrwise/Count";
		
		if(!new File(args[0]+"CBASE/Counts/").exists()){
			new File(args[0]+"CBASE/Counts/").mkdir();
		}
		if(!new File(args[0]+"CBASE/CountsChrwise/").exists()){
			new File(args[0]+"CBASE/CountsChrwise/").mkdir();
		}
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		try{
			FileInputStream in=new FileInputStream(file_samples);
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
			aa.add("PanCancer");
			entity_name=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entity_name[i]=aa.get(i);
			}
		}
		catch(Exception e){
			System.out.println(e);
		}		
		
		
		
		try{
			ArrayList<Gene>[] genes_chr=new ArrayList[chr.length];
			for (int i=0;i<genes_chr.length;i++){
				genes_chr[i]=new ArrayList<Gene>();
			}
			
			//read gene names with coordinates
			FileInputStream in=new FileInputStream(file_genes);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			String s="";
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
			
			ArrayList<Integer>[][] index_gene=new ArrayList[chr.length][];
			for (int i=0;i<chr.length;i++){
				index_gene[i]=new ArrayList[1+chr_length[i]/10000];
				for (int j=0;j<index_gene[i].length;j++){
					index_gene[i][j]=new ArrayList<Integer>();
				}
			}
			
			for (int i=0;i<genes_chr.length;i++){
				for (int j=0;j<genes_chr[i].size();j++){
					ArrayList<Integer> indices=new ArrayList<Integer>();
					for (int k=0;k<genes_chr[i].get(j).coord.size();k++){
						for (int p=genes_chr[i].get(j).coord.get(k)[0]/10000;p<=genes_chr[i].get(j).coord.get(k)[1]/10000;p++){
							if(!contains(p,indices)){
								indices.add(p);
							}
							
						}
					}
					for (int k=0;k<indices.size();k++){
						index_gene[i][indices.get(k)].add(j);
					}
				}
			}
			
			
			//initializie counters for the different categories
			 
			double[][][] lambda=new double[chr.length][][];//[genes.size()][entity_name.length];
			double[][][] lambda_syn=new double[chr.length][][];//[genes.size()][entity_name.length];
			int[][][] count=new int[chr.length][][];//[genes.size()][entity_name.length];
			int[][][] count_syn=new int[chr.length][][];//[genes.size()][entity_name.length];
			int[][][] count_maf=new int[chr.length][][];
			
			for (int i=0;i<chr.length;i++){
				lambda[i]=new double[genes_chr[i].size()][entity_name.length];
				lambda_syn[i]=new double[genes_chr[i].size()][entity_name.length];
			
				count[i]=new int[genes_chr[i].size()][entity_name.length];
				count_syn[i]=new int[genes_chr[i].size()][entity_name.length];
				count_maf[i]=new int[genes_chr[i].size()][entity_name.length];
				
			}
			
			//read the counts from the previous steps
			for (int j=0;j<entity_name.length;j++){
				//System.out.println(entity_name[j]);
				for (int i=0;i<chr.length;i++){
					in=new FileInputStream(file_count+entity_name[j]+"_Chr"+chr[i]+".txt");
					inn=new DataInputStream(in);
					input= new BufferedReader(new InputStreamReader(inn));
					while((s=input.readLine())!=null){
						String[] t=s.split("	");
						int ii=index_gene(t[0],genes_chr[i]);
						if(ii!=-1){
							lambda_syn[i][ii][j]=Double.parseDouble(t[4]);
							lambda[i][ii][j]=Double.parseDouble(t[3]);
							
							count_syn[i][ii][j]=Integer.parseInt(t[6]);
							count[i][ii][j]=Integer.parseInt(t[5]);
							
						}
					}
					input.close();
				}
			}
			
			
			Hashtable<String,Integer> table_entity=new Hashtable<String,Integer>();
			//ArrayList<String>[] list=new ArrayList[entity_name.length-1];
			//for (int i=0;i<list.length;i++){
			//	list[i]=new ArrayList<String>();
			//}	
			in=new FileInputStream(file_samples);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				table_entity.put(t[index_header[1]],index(t[index_header[2]],entity_name));
				
				//list[index(t[index_header[2]],entity_name)].add(t[index_header[1]]);
			}
			input.close();
			
			System.out.println("XXXX");
			//check counts from previous step against the counts in the maf file 
			//(eg. some counts in the intron exon boarder might not have been part
			//of the exon sequence but still be helpful to detected hypermutant genes
			in=new FileInputStream(file_maf);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header_m=index_header(input.readLine().split("	"),index_header_maf);
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				String ref=t[index_header_m[7]].toUpperCase();
				String tumor=t[index_header_m[9]].toUpperCase();
				if(tumor.equals(ref)||tumor.equals("")){
					tumor=t[index_header_m[8]].toUpperCase();
				}
				if(ref.length()!=1||tumor.length()!=1){
					continue;
				}
				if(!isNucleotide(ref)||!isNucleotide(tumor)){
					continue;
				}
				int chr_index=index(t[index_header_m[1]],chr);
				if(chr_index==-1){
					continue;
				}
				int sample_index=table_entity.get(t[index_header_m[10]]);//index(t[index_header_m[10]],list);
				ArrayList<Integer> gene_index=index(Integer.parseInt(t[index_header_m[2]]),genes_chr[chr_index],index_gene[chr_index][Integer.parseInt(t[index_header_m[2]])/10000]);
				for (int i=0;i<gene_index.size();i++){
					count_maf[chr_index][gene_index.get(i)][sample_index]++;
					count_maf[chr_index][gene_index.get(i)][entity_name.length-1]++;
				}
				
				
			}
			input.close();
			System.out.println("XXXX");
			
			//add counts that were missed in the exon sequence
			for (int i=0;i<count.length;i++){
				for (int j=0;j<count[i].length;j++){
					for (int k=0;k<count[i][j].length;k++){
						int diff=count_maf[i][j][k]-(count[i][j][k]+count_syn[i][j][k]);
						if(diff>0){
							count[i][j][k]+=diff;
						}
					}
				}
			}
			
			
			//output of counts separately for nonsyn and syn mutations in the CBASE format
			for (int k=0;k<entity_name.length;k++){
				
				FileWriter out1=new FileWriter(file_out1+entity_name[k]+".txt");
				BufferedWriter output1=new BufferedWriter(out1);
				output1.write("Name	TargetSize	Count");
				output1.newLine();
				for (int i=0;i<genes_chr.length;i++){
					for (int j=0;j<genes_chr[i].size();j++){
						if(lambda[i][j][k]*lambda_syn[i][j][k]>0){
							output1.write(genes_chr[i].get(j).name+"	"+lambda[i][j][k]+"	"+count[i][j][k]);//lambda.get(i)+"	"+count.get(i));
							output1.newLine();
						}
						
					}
				}
				output1.close();
				
				FileWriter out2=new FileWriter(file_out2+entity_name[k]+".txt");
				BufferedWriter output2=new BufferedWriter(out2);
				output2.write("Name	TargetSize	Count");
				output2.newLine();
				for (int i=0;i<genes_chr.length;i++){
					for (int j=0;j<genes_chr[i].size();j++){
						if(lambda[i][j][k]*lambda_syn[i][j][k]>0){
							output2.write(genes_chr[i].get(j).name+"	"+lambda_syn[i][j][k]+"	"+count_syn[i][j][k]);//lambda.get(i)+"	"+count.get(i));
							output2.newLine();
						}
						
					}
				}
				output2.close();
				
				for (int i=0;i<genes_chr.length;i++){
					FileWriter out3=new FileWriter(file_out3+entity_name[k]+"_Chr"+chr[i]+".txt");
					BufferedWriter output3=new BufferedWriter(out3);
					for (int j=0;j<genes_chr[i].size();j++){
						if(lambda[i][j][k]*lambda_syn[i][j][k]>0){
							output3.write(genes_chr[i].get(j).name+"	"+lambda[i][j][k]+"	"+lambda_syn[i][j][k]+"	"+count[i][j][k]+"	"+count_syn[i][j][k]);//lambda.get(i)+"	"+count.get(i));
							output3.newLine();
						}
						
					}
					output3.close();
				}
				
				
			}
			
		}
		catch(Exception e){
			System.out.println(e);
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
		}
	}
	
	public static boolean contains(int x, ArrayList<Integer> y){
		for (int i=0;i<y.size();i++){
			if(y.get(i)==x){
				return true;
			}
		}
		return false;
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
	
	public static ArrayList<Integer> index(int pos,ArrayList<Gene> genes, ArrayList<Integer> z){
		ArrayList<Integer> index=new ArrayList<Integer>();
		for(int ii=0;ii<z.size();ii++){
		//for (int i=0;i<genes.size();i++){
			int i=z.get(ii);
			if(genes.get(i).contains(pos)){
				index.add(i);
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
	
	public static boolean isNucleotide(String s){
		if(s.equals("A")){
			return true;
		}
		else if(s.equals("C")){
			return true;
		}
		else if(s.equals("G")){
			return true;
		}
		else if(s.equals("T")){
			return true;
		}
		return false;
	}
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	public static int index(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return i;
			}
		}
		return -1;
	}
	public static int index(String s, ArrayList<String>[] t){
		for (int i=0;i<t.length;i++){
			for (int j=0;j<t[i].size();j++){
				if(t[i].get(j).equals(s)){
					return i;
				}
			}
		}
		return -1;
	}
	
	private static class Gene{
		String name="";
		ArrayList<int[]> coord=new ArrayList<int[]>();
		//int start=-1;
		//int end=-1;
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
}
