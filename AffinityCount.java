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
 * 			sample its 20-nucleotide count vector, i.e.		*
 * 			it counts the occurence of each nucleotide		*
 * 			around its mutations. These count vectors are	*
 * 			needed for the subsequent clustering of the 	*
 * 			samples.										*
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

public class AffinityCount {
	static String file_annotation="";
	static String file_aligned="";
	static String file_samples="";
	static String file_out="";
	static String file_out2="";
	
	
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
		file_out=args[0]+"AffinityCounts/AffinityCount.txt";
		file_out2=args[0]+"AffinityCounts/TypeCount.txt";
		
		if(!new File(args[0]+"AffinityCounts/").exists()){
			new File(args[0]+"AffinityCounts/").mkdir();
		}
		
		
		try{
			Hashtable<String, Integer> sample_table=new Hashtable<String, Integer>();
			ArrayList<String> sample_list=new ArrayList<String>();
			//read sample table and link sample names with their index
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
			
			//calculate the count vector for each chr separately
			Count[] cc=new Count[chr.length];
			for (int i=0;i<chr.length;i++){
				cc[i]=count(chr[i]);
			}
			
			
			//sum up context count vector over chr
			int[][][][] count =new int[no_samples][20][6][4];
			for (int i=0;i<count.length;i++){
				for (int j=0;j<count[i].length;j++){
					for (int k=0;k<count[i][j].length;k++){
						for (int l=0;l<count[i][j][k].length;l++){
							for (int m=0;m<chr.length;m++){
								count[i][j][k][l]+=cc[m].count[i][j][k][l];
							}
						}
					}
				}
			}
			
			//sum up type count vector over chr
			int[][] count_type =new int[no_samples][6];
			for (int i=0;i<count_type.length;i++){
				for (int j=0;j<count_type[i].length;j++){
					for (int m=0;m<chr.length;m++){
						count_type[i][j]+=cc[m].count_type[i][j];
					}		
				}
			}
			
			//Output context count vector
			FileWriter out=new FileWriter(file_out);
			BufferedWriter output= new BufferedWriter(out);
			output.write("Sample");
			for (int j=0;j<count[0].length;j++){
				for (int k=0;k<count[0][j].length;k++){
					for (int l=0;l<count[0][j][k].length;l++){
						if(j<10){
							output.write("	Pos"+(j-10)+"_Type"+(k+1)+"_Nucl_"+new String[]{"A","C","G","T"}[l]);
						}
						else{
							output.write("	Pos"+(j-9)+"_Type"+(k+1)+"_Nucl_"+new String[]{"A","C","G","T"}[l]);
						}
					}
				}
			}
			output.newLine();
			
			
			for (int i=0;i<count.length;i++){
				output.write(sample_list.get(i));
				for (int j=0;j<count[i].length;j++){
					for (int k=0;k<count[i][j].length;k++){
						for (int l=0;l<count[i][j][k].length;l++){
							output.write("	"+count[i][j][k][l]);
						}
					}
				}
				output.newLine();
			}
			output.close();
			
			//Output Type count vector
			out=new FileWriter(file_out2);
			output= new BufferedWriter(out);
			output.write("Sample	Type1	Type2	Type3	Type4	Type5	Type6");
			output.newLine();
			for (int i=0;i<count_type.length;i++){
				output.write(sample_list.get(i));
				for (int j=0;j<count_type[i].length;j++){
					
					output.write("	"+count_type[i][j]);
					
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
	
	//compute the count vector on an individual chr
	public static Count count(String chr){
		int[][][][] count =new int[no_samples][20][6][4];
		int[][] count_type =new int[no_samples][6];
		try{
			ArrayList<Integer> queue_pos=new ArrayList<Integer>();
			ArrayList<String> queue_nucl=new ArrayList<String>();
			ArrayList<String> queue_sample1=new ArrayList<String>();
			ArrayList<String> queue_sample2=new ArrayList<String>();
			ArrayList<String> queue_sample3=new ArrayList<String>();
			
			String s="";
			
			//go through the aligned file from Step 1 and add positions into a queue
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
				
				//as soon as queue is large enoguh compute the nucleodide count around the center of the queue
				//update the count vectors and delete the first element of the queue
				if(queue_pos.size()>20){
					if(!queue_sample1.get(10).equals("")){
						String[] tt=queue_sample1.get(10).split(";");
						for(int i=-10;i<=10;i++){
							if(i==0){
								continue;
							}
							if(queue_pos.get(10)+i==queue_pos.get(10+i)){
								for (int j=0;j<tt.length;j++){
									if(queue_nucl.get(10).equals("C")||queue_nucl.get(10).equals("T")){
										count[Integer.parseInt(tt[j])][index_array(i)][type(queue_nucl.get(10),0)][index_nucl(queue_nucl.get(10+i))]++;
									}
									else{
										count[Integer.parseInt(tt[j])][19-index_array(i)][type(queue_nucl.get(10),0)][3-index_nucl(queue_nucl.get(10+i))]++;
									}
								}
								
							}
						}
						for (int j=0;j<tt.length;j++){
							count_type[Integer.parseInt(tt[j])][type(queue_nucl.get(10),0)]++;
						}
					}
					if(!queue_sample2.get(10).equals("")){
						String[] tt=queue_sample2.get(10).split(";");
						for(int i=-10;i<=10;i++){
							if(i==0){
								continue;
							}
							if(queue_pos.get(10)+i==queue_pos.get(10+i)){
								for (int j=0;j<tt.length;j++){
									if(queue_nucl.get(10).equals("C")||queue_nucl.get(10).equals("T")){
										count[Integer.parseInt(tt[j])][index_array(i)][type(queue_nucl.get(10),1)][index_nucl(queue_nucl.get(10+i))]++;
									}
									else{
										count[Integer.parseInt(tt[j])][19-index_array(i)][type(queue_nucl.get(10),1)][3-index_nucl(queue_nucl.get(10+i))]++;
									}
									
									//count[Integer.parseInt(tt[j])][index_array(i)][type(queue_nucl.get(10),1)][index_nucl(queue_nucl.get(10+i))]++;
								}
								
							}
						}
						for (int j=0;j<tt.length;j++){
							count_type[Integer.parseInt(tt[j])][type(queue_nucl.get(10),1)]++;
						}
					}
					if(!queue_sample3.get(10).equals("")){
						String[] tt=queue_sample3.get(10).split(";");
						for(int i=-10;i<=10;i++){
							if(i==0){
								continue;
							}
							if(queue_pos.get(10)+i==queue_pos.get(10+i)){
								for (int j=0;j<tt.length;j++){
									if(queue_nucl.get(10).equals("C")||queue_nucl.get(10).equals("T")){
										count[Integer.parseInt(tt[j])][index_array(i)][type(queue_nucl.get(10),2)][index_nucl(queue_nucl.get(10+i))]++;
									}
									else{
										count[Integer.parseInt(tt[j])][19-index_array(i)][type(queue_nucl.get(10),2)][3-index_nucl(queue_nucl.get(10+i))]++;
									}
									
									//count[Integer.parseInt(tt[j])][index_array(i)][type(queue_nucl.get(10),2)][index_nucl(queue_nucl.get(10+i))]++;
								}
								
							}
						}
						for (int j=0;j<tt.length;j++){
							count_type[Integer.parseInt(tt[j])][type(queue_nucl.get(10),2)]++;
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
		return new Count(count, count_type);
	}
	
	private static class Count{
		int[][][][] count=null;
		int[][] count_type=null;
		public Count(int[][][][] count, int[][] count_type){
			this.count=count;
			this.count_type=count_type;
		}
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
