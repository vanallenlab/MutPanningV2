/************************************************************           
 * MutPanning												*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This scripts starts from a Maf file and 		*
 * 			annotates each position in the human			*
 * 			exome with the sample indices which carry		*
 * 			a mutation. These aligned files are needed		*
 * 			for the subsequent steps of the algorithm		*										*   
 *************************************************************/


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;

public class AlignHG19 {

	static String file_annotation="";
	//static String file_peptide="";
	//static String file_coverage="";
	static String file_samples="";
	static String file_maf="";
	static String file_out="";
	//static String[] chr={};
	static String[] index_header_samples={"ID","Sample","Cohort"};
	static String[] index_header_maf={"Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode"};
	
	/*	argument 0: root file, where all the other files can be found
	 * argument 1: maf file
	 * argument 2: sample annotation file
	 * agrument 3: which chr should be analyzed
	 */
	
	
	//static Hashtable<Integer,Integer> position_table=new Hashtable<Integer,Integer>();
	//static String[] nucl=null;
	//static ArrayList<Position> position=new ArrayList<Position>();
	//static Position[] position=null;
	
	static String[] cc=new String[]{"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
	
	
	public static void main(String[] args){
		
		file_annotation=args[3]+"AnnotationHg19/Annotation_chr";
		//file_peptide=args[0]+"Hg19/ASAnnotationHg19/ASAnnotation_chr";
		//file_coverage=args[0]+"Hg19/CoverageExome_TCGA/Coverage_chr";
		file_samples=args[2];
		file_maf=args[1];
		file_out=args[0]+"AlignHg19/AlignHg19Chr";
		
		if(!new File(args[0]+"AlignHg19").exists()){
			new File(args[0]+"AlignHg19").mkdir();
		}
		
		/*String chr="";
		if(index(args[3],new String[]{"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"})!=-1){
			chr=args[3];
		}
		else{
			System.exit(0);
		}*/
		
	
		try{
			Hashtable<String, Integer> sample_table=new Hashtable<String, Integer>();
			// read sample names and link them to indices
			FileInputStream in=new FileInputStream(file_samples);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				sample_table.put(t[index_header[1]],Integer.parseInt(t[index_header[0]]));
			}
			input.close();
			
			int aa=0;
			in=new FileInputStream(file_maf);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				aa++;
			}
			input.close();
			
			String[] chrr=new String[aa];
			String[] ref=new String[aa];
			int[] pos=new int[aa];
			int[] index=new int[aa];
			int[] type=new int[aa];
			
			aa=0;
			in=new FileInputStream(file_maf);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header_m=index_header(input.readLine().split("	"),index_header_maf);
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				
				chrr[aa]=t[index_header_m[1]];
				pos[aa]=Integer.parseInt(t[index_header_m[2]]);
				ref[aa]=t[index_header_m[7]].toUpperCase();
				String tumor=t[index_header_m[9]].toUpperCase();
				if(tumor.equals(ref[aa])||tumor.equals("")){
					tumor=t[index_header_m[8]].toUpperCase();
				}
				
				index[aa]=sample_table.get(t[index_header_m[10]]);
				if(ref[aa].length()!=1||tumor.length()!=1){
					type[aa]=-1;
				}
				else if(!isNucleotide(ref[aa])||!isNucleotide(tumor)){
					type[aa]=-1;
				}
				else{
					type[aa]=type(ref[aa],tumor);
				}
				
				aa++;
				//if(aa%10000==0){
				//	System.out.println(aa);
				//}
			}
			input.close();
			
			for (int ll=0;ll<cc.length;ll++){
				System.out.println(cc[ll]);
				int nnn=0;
				in=new FileInputStream(file_annotation+cc[ll]+".txt");
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				while(input.readLine()!=null){
					nnn++;
				}
				input.close();
				
				String[] nucl=new String[nnn];
				Hashtable<Integer,Integer> position_table=new Hashtable<Integer,Integer>();
				
				in=new FileInputStream(file_annotation+cc[ll]+".txt");
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				String s1="";
				int kk=0;
				while((s1=input.readLine())!=null){
					String[] t1=s1.split("	");
					//if(kk%100000==0){
					//	System.out.println(kk+"/"+nnn);
					//}
					position_table.put(Integer.parseInt(t1[0]),kk);
					nucl[kk]=t1[1];
					kk++;
				}
				input.close();
				
				
				
				ArrayList<Integer> sample_index[]=new ArrayList[nucl.length];
				ArrayList<Integer> sample_type[]=new ArrayList[nucl.length];
				
				for (int k=0;k<pos.length;k++){
					
					if(!chrr[k].equals(cc[ll])){
						continue;
					}
					Integer ii=position_table.get(pos[k]);
					if(ii==null){
						continue;
					}
					
					if(type[k]==-1){
						continue;
					}
					
					if(!nucl[ii.intValue()].equals(ref[k])){
						continue;
					}
					
					if(sample_index[ii.intValue()]==null){
						sample_index[ii.intValue()]=new ArrayList<Integer>();
						sample_type[ii.intValue()]=new ArrayList<Integer>();
					}
					sample_index[ii.intValue()].add(index[k]);
					sample_type[ii.intValue()].add(type[k]);
				}
				
//				in=new FileInputStream(file_peptide+cc[ll]+".txt");
//				inn=new DataInputStream(in);
//				input= new BufferedReader(new InputStreamReader(inn));
//				
//				FileInputStream in2=new FileInputStream(file_coverage+cc[ll]+".txt");
//				DataInputStream inn2=new DataInputStream(in2);
//				BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
//				
				FileWriter out=new FileWriter(file_out+cc[ll]+".txt");
				BufferedWriter output= new BufferedWriter(out);
				
				s1="";
				String s2="";
				
				for (int j=0;j<nucl.length;j++){
					
					//s1=input.readLine();
					//s2=input2.readLine();
					//String[] t1=s1.split("	");
					//String[] t2=s2.split("	");
					
					//if(!t1[0].equals(t2[0])){
					//	System.exit(0);
					//}
					
					
					//output.write(t1[0]+"	"+t1[1]+"	"+t2[2]);
					int[] nn=new int[3];
					for (int k=0;k<3;k++){
						if(sample_index[j]!=null){
							for (int l=0;l<sample_index[j].size();l++){
								if(sample_type[j].get(l)==k){
									nn[k]++;
								}
								
							}
						}	
					}
					
					int i_max=-1;
					for (int i=0;i<nn.length;i++){
						if(nn[i]>0){
							i_max=i;
						}
					}
					
					if(i_max>=0){
						for (int k=0;k<=i_max;k++){
							if(k>0){
								output.write("	");
							}
							if(sample_index[j]!=null){
								int n=0;
								for (int l=0;l<sample_index[j].size();l++){
									if(sample_type[j].get(l)==k){
										if(n==0){
											output.write(""+sample_index[j].get(l));
										}
										else{
											output.write(";"+sample_index[j].get(l));
										}
										n++;
									}
									
								}
							}	
						}
					}
					
					//if(t1.length>=3){
					//	output.write("	"+t1[2]+"	"+t1[3]+"	"+t1[4]+"	"+t1[5]+"	"+t1[6]);
					//}
					output.newLine();
				
				}
				//input.close();
				//input2.close();
				
				
				output.close();	
				
				System.gc ();
				System.runFinalization ();
			}
			
	/*	
			int aa=0;
			ArrayList<Integer> sample_index[]=new ArrayList[nucl.length];
			ArrayList<Integer> sample_type[]=new ArrayList[nucl.length];
			//walk through the mutation annotation file. for each mutation link the index
			//of its sample to the genomic position of the reference genome
			
			in=new FileInputStream(file_maf);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header_m=index_header(input.readLine().split("	"),index_header_maf);
			while((s=input.readLine())!=null){
				//System.out.println(s);
				String[] t=s.split("	");
				aa++;
				if(aa%10000==0){
					System.out.println(aa);
				}
				if(!t[index_header_m[1]].equals(chr)){
					continue;
				}
				int pos=Integer.parseInt(t[index_header_m[2]]);
				Integer ii=position_table.get(pos);
				if(ii==null){
					continue;
				}
				
				
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
				
				int index=sample_table.get(t[index_header_m[10]]);
				int type=type(ref,tumor);
				if(type==-1){
					continue;
				}
				
				
				if(!nucl[ii.intValue()].equals(ref)){
					continue;
				}
				if(sample_index[ii.intValue()]==null){
					sample_index[ii.intValue()]=new ArrayList<Integer>();
					sample_type[ii.intValue()]=new ArrayList<Integer>();
				}
				sample_index[ii.intValue()].add(index);
				sample_type[ii.intValue()].add(type);
			}
			input.close();
			
		*/	
			/*
			in=new FileInputStream(file_maf);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header_m=index_header(input.readLine().split("	"),index_header_maf);
			while((s=input.readLine())!=null){
				//System.out.println(s);
				String[] t=s.split("	");
				aa++;
				if(aa%10000==0){
					System.out.println(aa);
				}
				if(!t[index_header_m[1]].equals(chr)){
					continue;
				}
				int pos=Integer.parseInt(t[index_header_m[2]]);
				Integer ii=position_table.get(pos);
				if(ii==null){
					continue;
				}
				
				
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
				
				int index=sample_table.get(t[index_header_m[10]]);
				int type=type(ref,tumor);
				if(type==-1){
					continue;
				}
				
				
				if(!nucl[ii.intValue()].equals(ref)){
					continue;
				}
				if(sample_index[ii.intValue()]==null){
					sample_index[ii.intValue()]=new ArrayList<Integer>();
					sample_type[ii.intValue()]=new ArrayList<Integer>();
				}
				sample_index[ii.intValue()].add(index);
				sample_type[ii.intValue()].add(type);
			}
			input.close();
			*/
			
			
			
			//now the position table contains all the positions together with the indices of the samples that contain mutations
			
			//output loop, generate a separate file for each chr
			
			
			
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
	
	public static int type(String ref, String tumor){//conversion to the type index
		
		if(ref.equals("A")){
			if(tumor.equals("C")){
				return 1;
			}
			else if(tumor.equals("G")){
				return 2;
			}
			else if(tumor.equals("T")){
				return 0;
			}
		}
		else if(ref.equals("C")){
			if(tumor.equals("A")){
				return 0;
			}
			else if(tumor.equals("G")){
				return 1;
			}
			else if(tumor.equals("T")){
				return 2;
			}
		}
		else if(ref.equals("G")){
			if(tumor.equals("A")){
				return 2;
			}
			else if(tumor.equals("C")){
				return 1;
			}
			else if(tumor.equals("T")){
				return 0;
			}
		}
		else if(ref.equals("T")){
			if(tumor.equals("A")){
				return 0;
			}
			else if(tumor.equals("C")){
				return 2;
			}
			else if(tumor.equals("G")){
				return 1;
			}
		}
		//System.out.println(ref+">"+tumor);
		return -1;
	}
	
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	
	//This method reads for 1 chromosome the reference sequence file
//	public static void readDELETE(String chr){
//		
//		try{
//			int nnn=0;
//			int nnn2=0;
//			
//			{
//				FileInputStream in=new FileInputStream(file_peptide+chr+".txt");
//				DataInputStream inn=new DataInputStream(in);
//				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
//				while(input.readLine()!=null){
//					nnn++;
//				}
//				input.close();
//			}
//			nucl=new String[nnn];
//			
//			FileInputStream in=new FileInputStream(file_peptide+chr+".txt");
//			DataInputStream inn=new DataInputStream(in);
//			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
//			String s1="";
//			int kk=0;
//			while((s1=input.readLine())!=null){
//				String[] t1=s1.split("	");
//				//if(kk%100000==0){
//				//	System.out.println(kk+"/"+nnn);
//				//}
//				position_table.put(Integer.parseInt(t1[0]),kk);
//				nucl[kk]=t1[1];
//				kk++;
//			}
//			input.close();
//			
//			/*
//			FileInputStream in=new FileInputStream(file_peptide+chr+".txt");
//			DataInputStream inn=new DataInputStream(in);
//			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
//			
//			FileInputStream in2=new FileInputStream(file_coverage+chr+".txt");
//			DataInputStream inn2=new DataInputStream(in2);
//			BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
//			
//			String s1="";
//			String s2="";
//			
//			int kk=0;
//			while((s1=input.readLine())!=null){
//				nnn2++;
//				if(nnn2%100000==0){
//					System.out.println(nnn2+"/"+nnn);
//				}
//				s2=input2.readLine();
//				String[] t1=s1.split("	");
//				String[] t2=s2.split("	");
//				
//				if(!t1[0].equals(t2[0])){
//					System.exit(0);
//				}
//				int pos=Integer.parseInt(t1[0]);
//				String nucl=t1[1];
//				double coverage=Double.parseDouble(t2[2]);
//				if(t1.length<3){
//					position_table.put(pos,kk);
//					//position.add(new Position(pos,nucl,coverage));
//					position[kk]=new Position(pos,nucl,coverage);
//				}
//				else{
//					position_table.put(pos,kk);
//					String as_ref=t1[2];
//					int as_no=Integer.parseInt(t1[3]);
//					String as_tumor1=t1[4];
//					String as_tumor2=t1[5];
//					String as_tumor3=t1[6];
//					//position.add(new Position(pos,nucl,coverage,as_ref,as_no,as_tumor1,as_tumor2,as_tumor3));
//					position[kk]=new Position(pos,nucl,coverage,as_ref,as_no,as_tumor1,as_tumor2,as_tumor3);
//				}
//				kk++;
//			}
//			input.close();
//			input2.close();
//			*/
//		}
//		catch(Exception e){
//			StackTraceElement[] aa=e.getStackTrace();
//			for (int i=0;i<aa.length;i++){
//				System.out.println(i+"	"+aa[i].getLineNumber());
//			}
//			System.out.println(e);
//		}
//	}
//	
	
	//the object encodes the annotation of each position and the indices of the samples which have a mutation
	private static class Position{
		int pos=-1;
		String nucl="";
		double coverage=0;
		String as_ref="";
		int as_no=-1;
		String as_tumor1="";
		String as_tumor2="";
		String as_tumor3="";
		ArrayList<Integer>[] samples=new ArrayList[3];
		
		
		public Position(){
			
		}
		
		public Position(int pos, String nucl, double coverage){
			this.pos=pos;
			this.nucl=nucl;
			this.coverage=coverage;
			samples[0]=new ArrayList<Integer>();
			samples[1]=new ArrayList<Integer>();
			samples[2]=new ArrayList<Integer>();
		}
		public Position(int pos, String nucl, double coverage,String as_ref, int as_no,String as_tumor1, String as_tumor2, String as_tumor3){
			this.pos=pos;
			this.nucl=nucl;
			this.coverage=coverage;
			this.as_ref=as_ref;
			this.as_no=as_no;
			this.as_tumor1=as_tumor1;
			this.as_tumor2=as_tumor2;
			this.as_tumor3=as_tumor3;
			samples[0]=new ArrayList<Integer>();
			samples[1]=new ArrayList<Integer>();
			samples[2]=new ArrayList<Integer>();
		}
		public void add(int sample_index,int type){
			samples[type].add(sample_index);
		}
	}
}
