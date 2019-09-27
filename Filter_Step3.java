/************************************************************           
 * MutPanning 												*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This is script is the final step to filter out	*
 * false positive significant genes. This steps writes all	*
 * genes which pass all filters into the final file of		*
 * significant genes. Further, this list transforms all		*
 * gene symbols into the HUGO nomenclature.					*
 * This is the final step of MutPanning						*
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
import java.util.Comparator;
import java.util.Hashtable;

public class Filter_Step3 {
	static String[] file_sign=new String[2];
	static String[] file_masked=new String[2];
	static String[] file_out=new String[2];
	
	
	
	static String[] entities=new String[0];
	static boolean[] compute_uniform=new boolean[0];
	static ArrayList<String> symbol=new ArrayList<String>();
	static Hashtable<String,Integer> table=new Hashtable<String,Integer>();
	static ArrayList<String> census=new ArrayList<String>(); 
	static ArrayList<String> black_list=new ArrayList<String>(); 
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	
	/*
	 * argument0: root file
	 * 
	 * argument1: file samples
	 */
	
	public static void main(String[] args, String[] args_entities, boolean[] args_compute_uniform){
		entities=args_entities;
		compute_uniform=args_compute_uniform;
		
		
		try{
			
			if(!new File(args[0]+"SignificanceFiltered/").exists()){
				new File(args[0]+"SignificanceFiltered/").mkdir();
			}
			file_sign=new String[]{args[0]+"SignificanceRaw/Significance",args[0]+"SignificanceRaw/SignificanceUniform"};
			file_masked=new String[]{args[0]+"PostSignFilter/MaskedPositions/MaskedPositions",args[0]+"PostSignFilter/MaskedPositions/MaskedPositionsUniform"};
			file_out=new String[]{args[0]+"SignificanceFiltered/Significance",args[0]+"SignificanceFiltered/SignificanceUniform"};
			
			//Read list of Approved Symbols. This list contains the HUGO symbol of each genes together with
			//all its aliases. Using this list all gene symbols are transformed to the HUGO nomenclature.
			FileInputStream in=new FileInputStream(args[2]+"SignificanceFilter/ApprovedSymbols.txt");
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			int n=0;
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				symbol.add(t[0]);
				table.put(t[0],n);
				String[] t1=t[3].split(", ");
				String[] t2=t[4].split(", ");
				if(t1.length>0&&!t1[0].equals("")){
					for (int i=0;i<t1.length;i++){
						table.put(t1[i],n);
					}
				}
				if(t2.length>0&&!t2[0].equals("")){
					for (int i=0;i<t2.length;i++){
						table.put(t2[i],n);
					}
				}
				n++;
			}
			input.close();
			
			//black listed genes are genes that produces false positive results in our pan cancer analysis
			in=new FileInputStream(args[2]+"SignificanceFilter/BlackList.txt");
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			while((s=input.readLine())!=null){
				black_list.add(symbol(s.split("	")[0]));
			}
			input.close();
			
			//read the cancer cenus as genes that should appear in the final report
			in=new FileInputStream(args[2]+"SignificanceFilter/census.csv");
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				census.add(symbol(s.split(",")[0]));
			}
			input.close();
			
		}
		catch(Exception e){
			System.out.println(e);
		}
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		/*
		try{
			FileInputStream in=new FileInputStream(args[1]);
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
			//aa.add("PanCancer");
			entities=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entities[i]=aa.get(i);
			}
		}
		catch(Exception e){
			System.out.println(e);
		}*/	
		
		String[] header_label={"Name","TargetSize","TargetSizeSyn","Count","CountSyn","SignificanceSyn","FDRSyn","Significance","FDR"};
		
		for (int l=0;l<2;l++){
			for (int k=0;k<entities.length;k++){
				try{
					if(!new File(file_sign[l]+entities[k]+".txt").exists()){
						continue;
					}
					
					//read the data from the raw significance gene file
					ArrayList<Gene> genes=new ArrayList<Gene>();
					FileInputStream in=new FileInputStream(file_sign[l]+entities[k]+".txt");
					DataInputStream inn=new DataInputStream(in);
					BufferedReader input= new BufferedReader(new InputStreamReader(inn));
					String[] header=input.readLine().split("	");
					int[] ii=new int[header_label.length];
					for (int i=0;i<header_label.length;i++){
						ii[i]=index(header_label[i],header);
					}
					
					
					String s="";
					while((s=input.readLine())!=null){
						String[] t=s.split("	");
						Gene gene=new Gene();
						gene.name=symbol(t[ii[0]]);
						gene.cov=Double.parseDouble(t[ii[1]]);
						gene.cov_syn=Double.parseDouble(t[ii[2]]);
						gene.count=Integer.parseInt(t[ii[3]]);
						gene.count_syn=Integer.parseInt(t[ii[4]]);
						gene.sign_complete_syn=Double.parseDouble(t[ii[5]]);
						gene.fdr_syn=Double.parseDouble(t[ii[6]]);
						gene.sign_complete=Double.parseDouble(t[ii[7]]);
						gene.fdr=Double.parseDouble(t[ii[8]]);
						
						//masking genes for which the synonymous (!!!) mutations diverge from the background distribution
						//if this occurs concordant deviation of the nonsynonymous mutations from the background distribution
						//does not necessarily indicate positive selection pressure, so that this criterion cannot be used
						//for cancer gene discovery in this case
						
						if(gene.fdr_syn<0.25){
							gene.pass=false;
						}
						genes.add(gene);
					}
					input.close();
					
					//masking all the genes that contain artifacts according to the BLAT filter
					if(new File(file_masked[l]+entities[k]+".txt").exists()){
						in=new FileInputStream(file_masked[l]+entities[k]+".txt");
						inn=new DataInputStream(in);
						input= new BufferedReader(new InputStreamReader(inn));
						while ((s=input.readLine())!=null){
							String[] t=s.split("	");
							genes.get(index_gene(symbol(t[0]),genes)).pass=false;
						}
						input.close();
					}
					
					
					for (int i=0;i<genes.size();i++){
						if(contains(genes.get(i).name,black_list)){
							genes.get(i).pass=false;
						}
						if(contains(genes.get(i).name,census)){
							genes.get(i).pass=true;
						}
					}
					
					
					
					Comparator<Gene> comp=(Gene g1, Gene g2)->{
						
						
						if(!g1.pass&&g2.pass){
							return +1;
						}
						if(g1.pass&&!g2.pass){
							return -1;
						}
						
						if(g1.fdr<g2.fdr){
							return -1;
						}
						if(g1.fdr>g2.fdr){
							return 1;
						}
						if(g1.count<g2.count){
							return +1;
						}
						if(g1.count>g2.count){
							return -1;
						}
						return 0;
					};
					
					
					Collections.sort(genes,comp);
					
					FileWriter out=new FileWriter(file_out[l]+entities[k]+".txt");
					BufferedWriter output= new BufferedWriter(out);
					//output.write("gene	target_n	target_s	count_n	count_s	p_seq	p_dm	p_cum	p	q");
					output.write("Name	TargetSize	TargetSizeSyn	Count	CountSyn	Significance	FDR");
					output.newLine();
					for (int i=0;i<genes.size();i++){
						if(genes.get(i).pass){//only output the genes that passed the filter
							output.write(genes.get(i).name+"	"+genes.get(i).cov+"	"+genes.get(i).cov_syn+"	"+genes.get(i).count+"	"+genes.get(i).count_syn+"	"+genes.get(i).sign_complete+"	"+genes.get(i).fdr);
							output.newLine();
						}
						
					}
					output.close();
					
				}
				catch(Exception e ){
					StackTraceElement[] aa=e.getStackTrace();
					for (int i=0;i<aa.length;i++){
						System.out.println(i+"	"+aa[i].getLineNumber());
					}
					System.out.println(e);
				}
			}
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
	
	public static String symbol(String s){
		int ii=index(s,symbol);
		if(ii!=-1){
			return s;
		}
		Integer a=table.get(s);
		if(a==null){
			return s;
		}
		
		return symbol.get(a);
		
	}
	
	public static int index_gene(String name, ArrayList<Gene> genes){
		for (int i=0;i<genes.size();i++){
			if(genes.get(i).name.equals(name)){
				return i;
			}
		}
		return -1;
	}
	
	private static class Gene{
		String name="";
		
		double cov=0;
		double cov_syn=0;
		int count=0;
		int count_syn=0;
		
		double sign_complete=1;
		double sign_complete_syn=1;
		
		double fdr=1;
		double fdr_syn=1;
		
		boolean pass=true;

		
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
