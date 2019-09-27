/************************************************************           
 * MutPanning 												*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This script determines for each individual		*
 * 			sample its 3-nucleotide count vector, i.e.		*
 * 			it counts the occurence of each nucleotide		*
 * 			around its mutations. These count vectors are	*
 * 			needed for the subsequent clustering of the 	*
 * 			samples. These count vectors are used for the	*
 * 			visualization of the mutaitonal processes		*
 * 			ongoing in each cluster to facilitate			*
 * 			cluster selection (COSMIC trinucleotide vector)	*
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
import java.util.Hashtable;

public class AffinityCount_Cosmic {
	static String file_annotation="";
	static String file_aligned="";
	static String file_samples="";
	static String file_out="";
	
	static int no_samples=0;
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	
	/*
	 * argument0: root file
	 * argument1: sample annotation file
	 */
	public static void main(String[] args){
		
		file_annotation=args[2]+"AnnotationHg19/Annotation_chr";
		file_aligned=args[0]+"AlignHg19/AlignHg19Chr";
		file_samples=args[1];
		file_out=args[0]+"AffinityCounts/CosmicCount.txt";
		
		if(!new File(args[0]+"AffinityCounts/").exists()){
			new File(args[0]+"AffinityCounts/").mkdir();
		}
		
		try{
			Hashtable<String, Integer> sample_table=new Hashtable<String, Integer>();
			ArrayList<String> sample_list=new ArrayList<String>();
			//read samples and link names to index
			FileInputStream in=new FileInputStream(file_samples);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				sample_table.put(t[index_header[1]],Integer.parseInt(t[index_header[0]]));
				sample_list.add(t[index_header[1]]);
			}
			input.close();
			no_samples=sample_list.size();
			//compute count vector on each chr separately
			
			
			//wait until all calculations are done
			int[][][] cc=new int[chr.length][][];
			for (int i=0;i<chr.length;i++){
				cc[i]=count(chr[i]);
			}
			
			//sum up count vectors over chr
			int[][] cosmic_count =new int[no_samples][96];
			for (int i=0;i<cosmic_count.length;i++){
				for (int j=0;j<cosmic_count[i].length;j++){
					for (int m=0;m<cc.length;m++){
						cosmic_count[i][j]+=cc[m][i][j];
					}		
				}
			}
			
			//output of trinucl count vectors
			FileWriter out=new FileWriter(file_out);
			BufferedWriter output= new BufferedWriter(out);
			output.write("Sample");
			for (int i=1;i<=96;i++){
				output.write("	Bin"+i);
			}
			output.newLine();
			for (int i=0;i<cosmic_count.length;i++){
				output.write(sample_list.get(i));
				for (int j=0;j<cosmic_count[i].length;j++){
					
					output.write("	"+cosmic_count[i][j]);
					
				}
				output.newLine();
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
	

	public static int[][] count(String chr){
		
		int[][]  cosmic_count =new int[no_samples][96];
		try{
			ArrayList<Integer> queue_pos=new ArrayList<Integer>();
			ArrayList<String> queue_nucl=new ArrayList<String>();
			ArrayList<String> queue_sample1=new ArrayList<String>();
			ArrayList<String> queue_sample2=new ArrayList<String>();
			ArrayList<String> queue_sample3=new ArrayList<String>();
			
			String s="";
			
			//go through the aligned mutation file and read positions into a queue
			FileInputStream in=new FileInputStream(file_aligned+chr+".txt");
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			
			FileInputStream in2=new FileInputStream(file_annotation+chr+".txt");
			DataInputStream inn2=new DataInputStream(in2);
			BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
		
			
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				String[] t2=input2.readLine().split("	");
				
				queue_pos.add(Integer.parseInt(t2[0]));
				queue_nucl.add(t2[1]);
				if(0<t.length){
					queue_sample1.add(t[0]);
				}
				else{
					queue_sample1.add("");
				}
				if(1<t.length){
					queue_sample2.add(t[1]);
				}
				else{
					queue_sample2.add("");
				}
				if(2<t.length){
					queue_sample3.add(t[2]);
				}
				else{
					queue_sample3.add("");
				}
				
				//as soon as queue is large engouth count nucleotides around the center of the queue
				//update count vectors and delete the first element in the queue
				if(queue_pos.size()>3){
					if(!queue_sample1.get(1).equals("")||!queue_sample2.get(1).equals("")||!queue_sample3.get(1).equals("")){
						if(queue_pos.get(1)+1==queue_pos.get(1+1)&&queue_pos.get(1)-1==queue_pos.get(1-1)){
							int triplet_index=-1;
							if(queue_nucl.get(1).equals("C")||queue_nucl.get(1).equals("T")){
								triplet_index=index_nucl(queue_nucl.get(0))*4+index_nucl(queue_nucl.get(2));
							}
							else if(queue_nucl.get(1).equals("G")||queue_nucl.get(1).equals("A")){
								triplet_index=(3-index_nucl(queue_nucl.get(2)))*4+(3-index_nucl(queue_nucl.get(0)));
							}
							
							String[] tt1=queue_sample1.get(1).split(";");
							String[] tt2=queue_sample2.get(1).split(";");
							String[] tt3=queue_sample3.get(1).split(";");
							
							if(queue_nucl.get(1).equals("C")||queue_nucl.get(1).equals("G")){
								if(!queue_sample1.get(1).equals("")){
									for (int j=0;j<tt1.length;j++){
										cosmic_count[Integer.parseInt(tt1[j])][0+triplet_index]++;
									}
								}
								if(!queue_sample2.get(1).equals("")){
									for (int j=0;j<tt2.length;j++){
										cosmic_count[Integer.parseInt(tt2[j])][16+triplet_index]++;
									}
								}
								if(!queue_sample3.get(1).equals("")){
									for (int j=0;j<tt3.length;j++){
										cosmic_count[Integer.parseInt(tt3[j])][32+triplet_index]++;
									}
								}
								
							}
							else if(queue_nucl.get(1).equals("A")||queue_nucl.get(1).equals("T")){
								if(!queue_sample1.get(1).equals("")){
									for (int j=0;j<tt1.length;j++){
										cosmic_count[Integer.parseInt(tt1[j])][48+0+triplet_index]++;
									}
								}
								if(!queue_sample2.get(1).equals("")){
									for (int j=0;j<tt2.length;j++){
										cosmic_count[Integer.parseInt(tt2[j])][48+32+triplet_index]++;//!!!!
									}
								}
								if(!queue_sample3.get(1).equals("")){
									for (int j=0;j<tt3.length;j++){
										cosmic_count[Integer.parseInt(tt3[j])][48+16+triplet_index]++;//!!!!
									}
								}
								
								
							}
						}
					
					}
					queue_pos.remove(0);
					queue_nucl.remove(0);
					queue_sample1.remove(0);
					queue_sample2.remove(0);
					queue_sample3.remove(0);
					
				}
			}
			input.close();
			input2.close();
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			System.out.println(e);
		}
		
		return cosmic_count;
	}
	
	public static int index_nucl(String s){
		if(s.equals("A")){
			return 0;
		}
		else if(s.equals("C")){
			return 1;
		}
		else if(s.equals("G")){
			return 2;
		}
		else if(s.equals("T")){
			return 3;
		}
		return -1;
	}
	
	public static int type(String n, int type){
		if(n.equals("C")||n.equals("G")){
			return type;
		}
		else if(n.equals("A")||n.equals("T")){
			if(type==0){
				return 3;
			}
			else if(type==1){
				return 5;
			}
			else if(type==2){
				return 4;
			}
		}
		return -1;
	}
	public static int index_array(int k){
		if(k<0){
			return k+10;
		}
		else if(k>0){
			return k+9;
		}
		return -1;
	}
	
}
