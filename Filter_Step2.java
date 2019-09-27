/************************************************************           
 * MutPanning 									*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This is script is the second step to filter out	*
 * false positive significant genes. This script evaluates	*
 * the output from the Blat queries of the previous step.	*
 * In brief, if the sequence context around a hotspot 		*
 * position was detected to have ambiguous possible			*
 * alignments this position and gene is masked from the 	*
 * significance test.										*
 *************************************************************/



import java.io.*;
import java.util.ArrayList;
import java.util.Collections;

public class Filter_Step2 {
	//static String[] file_query=new String[2];//"M:\\PostSignFilter\\Queries\\Query";
	//static String[] file_output=new String[2];//"M:\\PostSignFilter\\OutputsBLAT\\Output";
	//static String[] file_out=new String[2];//"M:\\PostSignFilter\\MaskedPositions";
	static String[] entities=new String[0];
	static boolean[] compute_uniform=new boolean[0];
	static String[] index_header_samples={"ID","Sample","Cohort"};
	static String[] types={"","Uniform"};
	static String file_out="";
	
	/*
	 * argument0: root file
	 * argument1: sample file 
	 */
	
	public static void main(String[] args, String[] args_entities, boolean[] args_compute_uniform){
		
		//file_query=new String[]{args[0]+"PostSignFilter/Queries/Query",args[0]+"PostSignFilter/Queries/QueryUniform"};
		//file_output=new String[]{args[0]+"PostSignFilter/OutputsBLAT/Output",args[0]+"PostSignFilter/OutputsBLAT/OutputUniform"};
		file_out=args[0]+"PostSignFilter/MaskedPositions/MaskedPositions";//new String[]{args[0]+"PostSignFilter/MaskedPositions/MaskedPositions",args[0]+"PostSignFilter/MaskedPositions/MaskedPositionsUniform"};
		String file_query=args[0]+"/PostSignFilter/Query.fa";
		String file_output=args[0]+"/PostSignFilter/OutputBLAT.txt";
		
		if(!new File(args[0]+"PostSignFilter/MaskedPositions/").exists()){
			new File(args[0]+"PostSignFilter/MaskedPositions/").mkdir();
		}
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		entities=args_entities;
		compute_uniform=args_compute_uniform;
		
		try{
			
			//for (int l=0;l<2;l++){
				//for (int k=0;k<entities.length;k++){
					/*if(!new File(file_query[l]+entities[k]+".fa").exists()){
						continue;
					}
					if(!new File(file_output[l]+entities[k]+".txt").exists()){
						continue;
					}*/
					
					
					//read all the queeries
					FileInputStream in=new FileInputStream(file_query);//[l]+entities[k]+".fa"
					DataInputStream inn=new DataInputStream(in);
					BufferedReader input= new BufferedReader(new InputStreamReader(inn));
					String s="";
					ArrayList<Query> query=new ArrayList<Query>();
					while((s=input.readLine())!=null){
						s=s.substring(1);//identifier of the query (including gene name)
						
						String seq=input.readLine();//sequence context of query
						query.add(new Query(s,seq));
					}
					input.close();
					
					int[] detected=new int[query.size()];
					
					
					//go through the query results
					//in brief associate the query results with the original queries
					//based on their identifier. test whether the alignment is identical
					//to the original sequence except for the hotspot mutation
					//if this is the case the hotspot mutation is likely to be a mis-
					//alignment artifact from a different region and is masked in the 
					//subsequent analysis
					in=new FileInputStream(file_output);//[l]+entities[k]+".txt"
					inn=new DataInputStream(in);
					input= new BufferedReader(new InputStreamReader(inn));
					while((s=input.readLine())!=null){
						if(s.contains("--------")){
							break;
						}
					}
					while((s=input.readLine())!=null){
						String[] t=s.split("	");
						if(!t[17].equals("1")){//no blocks
							continue;
						}
						String name_complete=t[9];
						//System.out.println(name_complete);
						int alignment_start=Integer.parseInt(t[11]);
						String original_seq=t[21].split(",")[0].toUpperCase();
						String aligned_seq=t[22].split(",")[0].toUpperCase();
						if(t[8].equals("-")){
							original_seq=reverse(original_seq);
							aligned_seq=reverse(aligned_seq);
						}
						//System.out.println(index(name_complete,query));
						if(index(name_complete,query)==-1){
							continue;
						}
						int[] alignment_index=new int[query.get(index(name_complete,query)).seq.length()];
						for (int i=0;i<alignment_index.length;i++){
							alignment_index[i]=-1;
						}
						for (int i=0;i<original_seq.length();i++){
							if(original_seq.charAt(i)==query.get(index(name_complete,query)).seq.charAt(alignment_start+i)){
								if(original_seq.charAt(i)==aligned_seq.charAt(i)){
									alignment_index[alignment_start+i]=1;
								}
								else{
									alignment_index[alignment_start+i]=0;
								}
							}
							else{
								alignment_index[alignment_start+i]=-1;
							}
						}
						
						String[] ttt=name_complete.split(":");
						String[] tt=ttt[1].split(";");
						String gene=tt[0].split("_")[0];
						String type=tt[0].split("_")[1];
						int count=Integer.parseInt(tt[2]);
						int start_index=Integer.parseInt(tt[1].split(",")[0]);
						int end_index=Integer.parseInt(tt[1].split(",")[1]);
						
						//System.out.println(name_complete+"	"+gene+"	"+type+"	"+count+"	"+start_index+"	"+end_index);
						
						if(!type.equals("alt")){
							continue;
						}
						if(alignment_index[-start_index]!=1){
							continue;
						}
						
						int prev=0;
						while(-start_index+prev>=0&&alignment_index[-start_index+prev]==1){
							prev--;
						}
						prev++;
						int post=0;
						while(-start_index+post<alignment_index.length&&alignment_index[-start_index+post]==1){
							post++;
						}
						post--;
						
						if(post-prev+1>=25){
							detected[index(name_complete,query)]=1;
						}
						
					}
					
					input.close();
					
					
					//for (int i=0;i<query.size();i++){
					//	System.out.println(query.get(i).name+"	"+detected[i]);
					//}
					
					FileWriter[][] out=new FileWriter[entities.length][2];//(file_out[l]+entities[k]+".txt");
					BufferedWriter[][] output= new BufferedWriter[entities.length][2];//(out);
					for (int i=0;i<entities.length;i++){
						out[i][0]=new FileWriter(file_out+entities[i]+".txt");
						output[i][0]=new BufferedWriter(out[i][0]);
						if(compute_uniform[i]){
							out[i][1]=new FileWriter(file_out+"Uniform"+entities[i]+".txt");
							output[i][1]=new BufferedWriter(out[i][1]);
						}
					}
					for (int i=0;i<query.size();i++){
						if(detected[i]==1){
							String[] ttt=query.get(i).name.split(":");
							String[] tttt=ttt[0].split("_");
							int kk=-1;
							if(tttt.length>1){
								kk=index(tttt[1],types);
							}
							else{
								kk=0;
							}
							//output[index(tttt[0],entities)][kk].write(query.get(i).name);
							//output[index(tttt[0],entities)][kk].newLine();
							output[index(tttt[0],entities)][kk].write(ttt[1].split(";")[0].split("_")[0]+"	"+ttt[1].split(";")[3].split(",")[0]+"	"+ttt[1].split(";")[3].split(",")[1]+"	"+ttt[1].split(";")[3].split(",")[2]);
							output[index(tttt[0],entities)][kk].newLine();
						}
					}
					for (int i=0;i<output.length;i++){
						for (int j=0;j<output[i].length;j++){
							if(output[i][j]!=null){
								output[i][j].close();
							}
						}
					}
					
				/*	
					//output of the misalginment positions found
					FileWriter out=new FileWriter(file_out[l]+entities[k]+".txt");
					BufferedWriter output= new BufferedWriter(out);
					for (int i=0;i<query.size();i++){
						if(detected[i]==1){
							output.write(query.get(i).name.split(";")[0].split("_")[0]+"	"+query.get(i).name.split(";")[3].split(",")[0]+"	"+query.get(i).name.split(";")[3].split(",")[1]+"	"+query.get(i).name.split(";")[3].split(",")[2]);
							output.newLine();
						}
					}
					output.close();
				*/
				//}
			//}
			
			
			
			
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			
			System.out.println(e);
		}
	}
	
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
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
	
	public static String ou(int[] array){
		String ou="";
		for (int i=0;i<array.length;i++){
			if(i>0){
				ou=ou+"|"+array[i];
			}
			else{
				ou=ou+array[i];
			}
			
		}
		return ou;
	}
	
	public static String reverse(String s){
		String rev="";
		for (int i=0;i<s.length();i++){
			if(s.charAt(i)=='A'){
				rev="T"+rev;
			}
			else if(s.charAt(i)=='C'){
				rev="G"+rev;
			}
			else if(s.charAt(i)=='G'){
				rev="C"+rev;
			}
			else if(s.charAt(i)=='T'){
				rev="A"+rev;
			}
			else{
				rev="N"+rev;
			}
		}
		return rev;
	}
	
	public static int index(String name, ArrayList<Query> query){
		for (int i=0;i<query.size();i++){
			if(query.get(i).name.equals(name)){
				return i;
			}
		}
		return -1;
	}
	
	private static class Query{
		String name="";
		String seq="";
		public Query(String name, String seq){
			this.name=name;
			this.seq=seq;
		}
		
	}
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
}
