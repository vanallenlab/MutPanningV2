/************************************************************           
 * MutPanning - Step 6									*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*   
 *															*   
 * Summary: This script first computes the mutation rate of	*
 * 			each position in the human exome, based on the	*
 * 			clusters that were dervied in the previous 		*
 * 			steps. Based on these mutation rates, the script*
 * 			computes the mutational significance of each 	*
 * 			gene. The major component of this test is based *
 * 			on the combined null hypothesis model, which	*
 * 			integrates mutation counts and sequence context	*
 *  		around mutations. Further this method tests for	*
 *  		local clustering of mutations in mutation 		*
 *  		hotspots (enhancing the sensitivity for 		*
 *  		oncogenes)and for accumulation of protein 		*
 *  		damaging mutations (enhancing the sensitivity 	*
 *  		for tumor supressor genes). With the latter test*
 *  		insertions and deletions are considered in our 	*
 *  		test statistics.								*
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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.SplittableRandom;
import java.util.concurrent.ThreadLocalRandom;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.special.Gamma;
import jdistlib.math.Bessel;





public class ComputeSignificance {
	static ArrayList<double[][][]> lambda_context=new ArrayList<double[][][]>();
	static ArrayList<double[][]> lambda_type=new ArrayList<double[][]>();
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
	
	static double[][][] lambda_context_product6_weight=null;
	
	static ArrayList<Integer> pos=null;
	static ArrayList<String> nucl=null;
	static ArrayList<Integer> nucl_index=null;
	static ArrayList<Double> coverage=null;
	static ArrayList<int[]> label=null;
	static ArrayList<int[][]> count=null;
	static ArrayList<String> amino_acid=null;
	
	
	
	static String[] entities=null;
	static String[] index_header_samples={"ID","Sample","Cohort"};
	static int no_clusters=-1;
	static ArrayList<Gene>[] genes=new ArrayList[chr.length];
	//static ArrayList<Integer>[] list=new ArrayList[0];
	static Hashtable<Integer,Integer> table_entity=null;
	
	static String file_out="";
	static String file_out_uniform="";
	static String file_annotation="";
	static String file_align="";
	static String file_signatures="";
	static String file_reference="";
	static String file_clusters="";
	static String file_type="";
	static String file_samples="";
	static String file_destructive="";
	static String file_genes="";
	//static String file_out2="";
	static String file_count_genes="";
	static String file_table_bessel="";
	
	static double[][] weights=null;
	static ArrayList<Integer>[] weights_index=null;
	static double[][] params=null;
	static int[] mod_C=null;
	static double[] gamma_precompute=new double[10000];
	static String[] numbers=new String[]{"0","1","2","3","4","5","6","7","8","9"};
	
	//static int ii=-1;
	/*
	 * argument0: root file
	 * 
	 */
	//static int iii=-1;
	
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
			"ASMT",//double genes
			"ASMTL",
			"CSF2RA",
			"DHRSX",
			"IL3RA",
			"IL3R",
			"IL9R",
			"SPRY3",
			"ZBED1"
		};
	static Hashtable<String,Integer> table_nucl=new Hashtable<String,Integer>();
	
	static ThreadLocalRandom random_conc=null;
	//static SplittableRandom random_conc=null;
	
	
	//static Hashtable<String,Integer> table_failed_genes=new  Hashtable<String,Integer>();
	static double[][] table_bessel=new double[10000][250];
	static int[] non_syn_count=null;
	static double non_syn_coverage=0;
	static boolean[] compute_uniform=null;
	
	
	public static void main(String[] args, String[] args_entities, boolean[] args_compute_uniform){
//		for (int k=0;k<=100;k++){
//			System.out.println(k+"	"+probability_binom(100,k,0.005));
//		}
//		
//		System.exit(0);
		
		compute_uniform=args_compute_uniform;
		entities=args_entities;
		
		//System.exit(0);
		
		/*
		{
			//double x1=20;
			double x2=8.08;
			for (double n=1;n<=500;n++){
				
				System.out.println(n+"	"+Math.log(Bessel.k(x2, n, false))+"	"+Math.log(Bessel.k(x2, n, true)));
				//+"	"+Bessel.k(x1, n, false)
			}
			
			for (double alpha=1;alpha<=1000;alpha*=1.1){
				System.out.println(alpha+"	"+Math.log(Bessel.k(alpha, 300, false))+"	"+Math.log(Bessel.k(alpha, 300, true)));
			}
			
			System.exit(0);
		}
		*/
		
		/*
		for (int i=0;i<failed_genes.length;i++){
			table_failed_genes.put(failed_genes[i],i);
		}
		*/
		
		random_conc=ThreadLocalRandom.current();//new SplittableRandom();//
		file_out=args[0]+"SignificanceRaw/Significance";
		file_out_uniform=args[0]+"SignificanceRaw/SignificanceUniform";
		
		//file_out=args[0]+"MutationRateClusters/Lambda_Chr";
		file_annotation=args[2]+"AnnotationHg19/Annotation_chr";
		file_align=args[0]+"AlignHg19/AlignHg19Chr";
		file_signatures=args[0]+"ClusteringComplete/ClusteringComplete_Affinity.txt";
		file_reference=args[2]+"FileReferenceCount.txt";
		file_table_bessel=args[2]+"TableBessel.txt";
		file_clusters=args[0]+"ClusteringComplete/ClusteringComplete_Samples.txt";
		file_type=args[0]+"AffinityCounts/TypeCount.txt";
		file_destructive=args[0]+"CountDestructive/";
		file_samples=args[1];
		file_genes=args[2]+"Exons_Hg19.txt";
		//file_out2=args[0]+"CBASE/CountsRaw/Count";
		file_count_genes=args[0]+"CBASE/CountsChrwise/Count";
		
		
		if(!new File(args[0]+"SignificanceRaw/").exists()){
			new File(args[0]+"SignificanceRaw/").mkdirs();
		}
		
		//if (!new File(args[0]+"CBASE/CountsRaw/").exists()){
		//	new File(args[0]+"CBASE/CountsRaw/").mkdirs();
		//}
		
		table_nucl.put("A",0);
		table_nucl.put("C",1);
		table_nucl.put("G",2);
		table_nucl.put("T",3);
		table_nucl.put("a",0);
		table_nucl.put("c",1);
		table_nucl.put("g",2);
		table_nucl.put("t",3);
		
		try{
			//if(!new File(args[0]+"MutationRateClusters/").exists()){
			//	new File(args[0]+"MutationRateClusters/").mkdir();
			//}
			
			FileInputStream in=new FileInputStream(file_table_bessel);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			for (int i=0;i<table_bessel.length;i++){
				String[] t=input.readLine().split("	");
				for (int j=0;j<table_bessel[i].length;j++){
					table_bessel[i][j]=Double.parseDouble(t[j+1]);
				}
			}
			input.close();
			
			//each cluster delivers a "signature" for the distribution of passenger
			//mutations in that cluster. read that signature
			ArrayList<double[][][]> frequency=new ArrayList<double[][][]>();
			in=new FileInputStream(file_signatures);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int[][][] a =new int [20][6][4];
				int n=0;
				for (int i=0;i<a.length;i++){
					for (int j=0;j<a[i].length;j++){
						for (int k=0;k<a[i][j].length;k++){
							a[i][j][k]=Integer.parseInt(t[n+1]);
							n++;
						}
					}
				}
				int[] b=new int[6];
				int sum=0;
				for (int i=0;i<6;i++){
					for (int j=0;j<4;j++){
						b[i]+=a[9][i][j]+a[10][i][j];
						sum+=a[9][i][j]+a[10][i][j];
						//System.out.println(sum);
					}
				}
				frequency.add(freq(a));
				double[][] lambda_t=new double[2][];
				lambda_t[0]=new double[]{2*(double)(b[0]+b[1]+b[2])/(double)(sum),2*(double)(b[3]+b[4]+b[5])/(double)(sum)};
				lambda_t[1]=new double[]{3*(double)(b[0]+b[3])/(double)(sum),3*(double)(b[1]+b[5])/(double)(sum),3*(double)(b[2]+b[4])/(double)(sum)};
				lambda_type.add(lambda_t);
			}
			input.close();
			
			
			for (int i=0;i<gamma_precompute.length;i++){
				gamma_precompute[i]=Gamma.logGamma(i);
			}
			
			//read the distribution in the reference genome
			int[][][] reference_abs=new int[20][2][4];
			in=new FileInputStream(file_reference);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			for (int i=0;i<2;i++){
				for (int j=0;j<4;j++){
					String[] t=input.readLine().split("	");
					for (int k=0;k<20;k++){
						if(k<10){
							reference_abs[k][i][j]=Integer.parseInt(t[k+2]);
						}
						else{
							reference_abs[k][i][j]=Integer.parseInt(t[k+3]);
						}
					}
				}
			}
			input.close();
			
			//Compute the likelihood ratios lambda as the ratio between reference and observed counts
			double[][][] reference=freq(reference_abs);
			for (int i=0;i<frequency.size();i++){
				double[][][] ll=new double[20][6][4];
				for (int j=0;j<frequency.get(i).length;j++){
					for (int k=0;k<frequency.get(i)[j].length;k++){
						for (int l=0;l<frequency.get(i)[j][k].length;l++){
							ll[j][k][l]=frequency.get(i)[j][k][l]/reference[j][k/3][l];
						}
					}
				}
				lambda_context.add(ll);
			}
			
			int[] ttt={0,1,2,0,2,1};
			
//			lambda_context_product1=new double[lambda_context.size()][6][4][4][4][4][4];
//			lambda_context_product2=new double[lambda_context.size()][6][4][4][4][4][4];
//			lambda_context_product3=new double[lambda_context.size()][6][4][4][4][4][4];
//			lambda_context_product4=new double[lambda_context.size()][6][4][4][4][4][4];
			
			
			
			
			
			/*
			double[][][] lambda_context_product6=new double[lambda_context.size()][6][(int)(Math.pow(4,10))];
				//System.out.println(a+"/"+lambda_context.size());
				for (int k=0;k<6;k++){
					int[] x=new int[10];
					for (x[0]=0;x[0]<4;x[0]++){
						//System.out.println(x[0]);
						for (x[1]=0;x[1]<4;x[1]++){
							for (x[2]=0;x[2]<4;x[2]++){
								for (x[3]=0;x[3]<4;x[3]++){
									for (x[4]=0;x[4]<4;x[4]++){
										for (x[5]=0;x[5]<4;x[5]++){
											for (x[6]=0;x[6]<4;x[6]++){
												for (x[7]=0;x[7]<4;x[7]++){
													for (x[8]=0;x[8]<4;x[8]++){
														for (x[9]=0;x[9]<4;x[9]++){
															int len=0;
															int n=1;
															for (int l=0;l<x.length;l++){
																len+=x[l]*n;
																n*=4;
															}
															for (int a=0;a<lambda_context.size();a++){
																double prod=lambda_type.get(a)[0][k/3]*lambda_type.get(a)[1][ttt[k]];
																for (int l=0;l<x.length;l++){
																	prod*=lambda_context.get(a)[l+5][k][x[l]];	
																}
																lambda_context_product6[a][k][len]=prod;
																System.out.println(lambda_context_product6[a][k][len]+"	"+lambda_context_product6X[a][k][len]);
															}
														
														}
													}
												}
											}
										}
										
								
									}
								}
							}
						}
					}
				}
				
			
			*/
			
			
			System.out.println(System.currentTimeMillis());
			
			/*
			for (int a=0;a<lambda_context.size();a++){
				
				for (int k=0;k<6;k++){
					int[] x=new int[5];
					for (x[0]=0;x[0]<4;x[0]++){
						for (x[1]=0;x[1]<4;x[1]++){
							for (x[2]=0;x[2]<4;x[2]++){
								for (x[3]=0;x[3]<4;x[3]++){
									for (x[4]=0;x[4]<4;x[4]++){
										
										
										
										lambda_context_product1[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]=lambda_type.get(a)[0][k/3]*lambda_type.get(a)[1][ttt[k]];
										lambda_context_product2[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]=1;
										lambda_context_product3[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]=1;
										lambda_context_product4[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]=1;
										
										for (int l=0;l<x.length;l++){
											lambda_context_product1[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]*=lambda_context.get(a)[l][k][x[l]];
											lambda_context_product2[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]*=lambda_context.get(a)[l+5][k][x[l]];
											lambda_context_product3[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]*=lambda_context.get(a)[l+10][k][x[l]];
											lambda_context_product4[a][k][x[0]][x[1]][x[2]][x[3]][x[4]]*=lambda_context.get(a)[l+15][k][x[l]];
										}

										
									}
								}
							}
						}
					}
				}
				
			}*/
			
			
			/*
			in=new FileInputStream(file_samples);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
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
			}*/
			
			ArrayList<String> names_count=new ArrayList<String>(); 
			ArrayList<Integer> count=new ArrayList<Integer>();
			
			//read the no. mutaions per samples
			in=new FileInputStream(file_type);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				names_count.add(t[0]);
				int c=0;
				for (int i=1;i<=6;i++){
					c+=Integer.parseInt(t[i]);
				}
				count.add(c);
			}
			input.close();
			
			no_clusters=0;
			in=new FileInputStream(file_clusters);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				if(Integer.parseInt(s.split("	")[2])>no_clusters){
					no_clusters=Integer.parseInt(s.split("	")[2]);
				}
			}
			no_clusters++;
			input.close();
			
			weights=new double[entities.length+1][no_clusters];
			
			//sum up the no. mutations in each cluster
			System.out.println("START");
			in=new FileInputStream(file_clusters);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			input.readLine();
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int c=count.get(index(t[0],names_count));
				int index_entity=index(t[1],entities);
				int cluster=Integer.parseInt(t[2]);
				
				weights[index_entity][cluster]+=c;
				weights[weights.length-1][cluster]+=c;
			}
			input.close();
			
			
			//normalize weights. weight <0.01 are set to 0
			//this safes a lot of run time and does not affect the
			//mutation rate of the sample substantially
			for (int i=0;i<weights.length;i++){
				double sum=0;
				for (int j=0;j<weights[i].length;j++){
					sum+=weights[i][j];
				}
				for (int j=0;j<weights[i].length;j++){
					weights[i][j]/=sum;
				}
				for (int j=0;j<weights[i].length;j++){
					if(weights[i][j]<0.01){
						weights[i][j]=0;
					}
				}
				sum=0;
				for (int j=0;j<weights[i].length;j++){
					sum+=weights[i][j];
				}
				for (int j=0;j<weights[i].length;j++){
					weights[i][j]/=sum;
				}
			}
			
			weights_index=new ArrayList[weights.length];
			for (int i=0;i<weights.length;i++){
				weights_index[i]=new ArrayList<Integer>();
				for (int j=0;j<weights[i].length;j++){
					if(weights[i][j]>0){
						weights_index[i].add(j);
					}
				}
			}
			
			System.out.println("Step2	"+System.currentTimeMillis());
			lambda_context_product6_weight=new double[weights.length][6][(int)(Math.pow(4, 10))];
			
			for (int k=0;k<6;k++){
				double[][] lambda_context_product6=new double[lambda_context.size()][1];
				
				for (int a=0;a<lambda_context.size();a++){
					lambda_context_product6[a][0]=lambda_type.get(a)[0][k/3]*lambda_type.get(a)[1][ttt[k]];
					for (int l=0;l<10;l++){
						double[] array_new=new double[4*lambda_context_product6[a].length];
						for (int m=0;m<array_new.length;m++){
							array_new[m]=lambda_context_product6[a][m/4]*lambda_context.get(a)[14-l][k][m%4];//l+5
						}
						lambda_context_product6[a]=array_new;
						System.gc();
						System.runFinalization();
					}	
				}
				
				for (int i=0;i<lambda_context_product6_weight.length;i++){
					for (int j=0;j<lambda_context_product6[0].length;j++){
						for (int a=0;a<weights_index[i].size();a++){
							lambda_context_product6_weight[i][k][j]+=weights[i][weights_index[i].get(a)]*lambda_context_product6[weights_index[i].get(a)][j];
						}		
					}
					
				}
			}
			
			System.out.println(System.currentTimeMillis());
			
			for (int i=0;i<lambda_context_product6_weight.length;i++){
				for (int j=0;j<lambda_context_product6_weight[i].length;j++){
					for (int k=0;k<lambda_context_product6_weight[i][j].length;k++){
						if(lambda_context_product6_weight[i][j][k]>1000){
							lambda_context_product6_weight[i][j][k]=1000;
						}
					}
				}
			}
			
			System.out.println("Step2	"+System.currentTimeMillis());
			
			
			
			
			for (int i=0;i<genes.length;i++){
				genes[i]=new ArrayList<Gene>();
			}
				
			//read genes with positions
			in=new FileInputStream(file_genes);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				if(index(t[0],exclude)!=-1){
					continue;
				}
				int ii=index_gene(t[0],genes[Integer.parseInt(t[1])-1]);
				if(ii==-1){
					genes[Integer.parseInt(t[1])-1].add(new Gene(t[0],Integer.parseInt(t[2]),Integer.parseInt(t[3]),entities.length));
				}
				else{
					genes[Integer.parseInt(t[1])-1].get(ii).coord.add(new int[]{Integer.parseInt(t[2]),Integer.parseInt(t[3])});
				}
			}
			input.close();
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					genes[i].get(j).start=min(genes[i].get(j).coord);
					genes[i].get(j).end=max(genes[i].get(j).coord);
				}
			}
			
			
			/*for (int i=0;i<genes.length;i++){
				for (int j=genes[i].size()-1;j>=0;j--){
					if(genes[i].get(j).name.length()>2&&genes[i].get(j).name.substring(0, 2).equals("OR")&&contains(genes[i].get(j).name.substring(2,3),numbers)){//&&genes[i].get(j).ell_s>0.0
						genes[i].remove(j);
					}
				}
			}*/
			//for (int i=1;i<genes.length;i++){
			//	genes[i]=new ArrayList<Gene>();
			//}
			
			
			Hashtable<String,int[]> table_gene=new Hashtable<String,int[]>();
			Hashtable<String,Integer>[] table_gene2=new Hashtable[chr.length];
			for (int i=0;i<chr.length;i++){
				table_gene2[i]=new Hashtable<String,Integer>();
			}
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					table_gene.put(genes[i].get(j).name,new int[]{i,j});
					table_gene2[i].put(genes[i].get(j).name,j);
				}
			}
			
			
			table_entity=new Hashtable<Integer,Integer>();
			//list=new ArrayList[entities.length];
			//for (int i=0;i<list.length;i++){
			//	list[i]=new ArrayList<Integer>();
			//}
			
			in=new FileInputStream(file_samples);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int ii=index(t[index_header[2]],entities);
				table_entity.put(Integer.parseInt(t[index_header[0]]),ii);
				//list[ii].add(Integer.parseInt(t[index_header[0]]));
			}
			input.close();
			
			mod_C=new int[entities.length];
			params=new double[entities.length][];
			for(int k=0;k<entities.length;k++){
				for (int i=0;i<chr.length;i++){
					in=new FileInputStream(file_count_genes+entities[k]+"_Chr"+chr[i]+".txt");
					inn=new DataInputStream(in);
					input= new BufferedReader(new InputStreamReader(inn));
					
					while((s=input.readLine())!=null){
						String[] t=s.split("	");
						int ii=index_gene(t[0],genes[i]);
						if(ii!=-1){
							genes[i].get(ii).cov=Double.parseDouble(t[1]);
							genes[i].get(ii).cov_syn=Double.parseDouble(t[2]);
							genes[i].get(ii).count[k]=Integer.parseInt(t[3]);
							genes[i].get(ii).count_syn[k]=Integer.parseInt(t[4]);
						}
						
					}
					input.close();
				}
				
				
				String infile_syn = args[0]+"CBASE/Counts/CountSilent"+entities[k]+".txt";		//	synonymous somatic mutation data input file
				String infile_nonsyn = args[0]+"CBASE/Counts/Count"+entities[k]+".txt";		//	nonsynonymous somatic mutation data input file
				String infile_param = args[0]+"CBASE/Parameters_Summary/Parameters"+entities[k]+".txt";
				
				in=new FileInputStream(infile_syn);
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				input.readLine();
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					int[] ii=table_gene.get(t[0]);
					if(ii!=null){
						genes[ii[0]].get(ii[1]).ell_s[k]=Double.parseDouble(t[1]);
						genes[ii[0]].get(ii[1]).sobs[k]=Double.parseDouble(t[2]);
					}
				}
				input.close();
				
				in=new FileInputStream(infile_nonsyn);
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				input.readLine();
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					int[] ii=table_gene.get(t[0]);
					if(ii!=null){
						genes[ii[0]].get(ii[1]).ell_x[k]=Double.parseDouble(t[1]);
						genes[ii[0]].get(ii[1]).xobs[k]=Double.parseDouble(t[2]);
					}
				}
				input.close();
				
				
				
				
				ArrayList<double[]> all_models=new ArrayList<double[]>();
				in=new FileInputStream(infile_param);//
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				while((s=input.readLine())!=null){
					String[] field=s.split(", ");
					double[] f=new double[field.length];
					for (int i=0;i<field.length;i++){
						f[i]=Double.parseDouble(field[i]);
					}
					all_models.add(f);
				}
				input.close();
				
				int[] modC_map = {2,2,4,4,5,5};						// map model --> number of params

				double cur_min=1e20; 
				int cur_ind=10;
				for (int m=0;m<all_models.size();m++){
					if (2.*modC_map[(int)(all_models.get(m)[all_models.get(m).length-1])-1] + 2.*all_models.get(m)[all_models.get(m).length-2] < cur_min){
						cur_min = 2.*modC_map[(int)(all_models.get(m)[all_models.get(m).length-1])-1] + 2.*all_models.get(m)[all_models.get(m).length-2];
						cur_ind = m;
					}
				}

				mod_C[k]	= (int)(all_models.get(cur_ind)[all_models.get(cur_ind).length-1]);
				params[k]	= sub(all_models.get(cur_ind),0,all_models.get(cur_ind).length-2);//[:-2];
				
				//System.out.println("Using model "+mod_choice[mod_C]);
//				System.out.print(entities[k]+"	"+mod_C[k]);
//				for (int i=0;i<params[k].length;i++){
//					System.out.print("	"+params[k][i]);
//				}
//				System.out.println();
			}
			//System.exit(0);
			
			/*
			for (int i=0;i<weights.length;i++){
				for (int j=0;j<weights[i].length;j++){
					System.out.print("	"+weights[i][j]);
				}
				System.out.println();
			}
			System.exit(0);*/
			
			//for (int i=0;i<chr.length;i++){
			//	run(i);
			//}
			
			/*
			for (int i=0;i<lambda_context.size();i++){
				for (int j=0;j<lambda_context.get(i).length;j++){
					for (int k=0;k<lambda_context.get(i)[j].length;k++){
						for (int l=0;l<lambda_context.get(i)[j][k].length;l++){
							System.out.print("	"+lambda_context.get(i)[j][k][l]);
						}
						
					}
				}
				System.out.println();
			}
			System.exit(0);
			*/
			//iii=index("Skin",entities);
			
			
			non_syn_count=new int[entities.length];
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					for (int k=0;k<entities.length;k++){
						non_syn_count[k]+=genes[i].get(j).count[k];
					}
					non_syn_coverage+=genes[i].get(j).cov;
				}
			}
			
			
			
			
			System.out.println("START");
			for (int i=0;i<chr.length;i++){
				//if(chr[i].equals("3")){
				run(i);
				//}
				//run(index("Y",chr));
			}
			
			/*
			for (int i=0;i<entities.length;i++){
				in=new FileInputStream(file_out+entities[i]+"V2.txt");
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				input.readLine();
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					int[] ii=table_gene.get(t[0]);
					if(ii!=null){
						genes[ii[0]].get(ii[1]).sign_vector_syn[i]=Double.parseDouble(t[5]);
						genes[ii[0]].get(ii[1]).sign_hotspot_syn[i]=1;//Double.parseDouble(t[6]);
						//genes[ii[0]].get(ii[1]).sign_complete_syn[i]=Double.parseDouble(t[7]);ddd
						genes[ii[0]].get(ii[1]).sign_combined[i]=Double.parseDouble(t[8]);
						//genes[ii[0]].get(ii[1]).sign_destruct[i]=Double.parseDouble(t[9]);
						genes[ii[0]].get(ii[1]).sign_hotspot[i]=Double.parseDouble(t[10]);
						//genes[ii[0]].get(ii[1]).sign_complete[i]=Double.parseDouble(t[11]);
					}
					
				}
				input.close();
				
				if(compute_uniform[i]){
					in=new FileInputStream(file_out_uniform+entities[i]+"V2.txt");
					inn=new DataInputStream(in);
					input= new BufferedReader(new InputStreamReader(inn));
					input.readLine();
					while((s=input.readLine())!=null){
						String[] t=s.split("	");
						int[] ii=table_gene.get(t[0]);
						if(ii!=null){
							//genes[ii[0]].get(ii[1]).sign_vector_syn[i]=Double.parseDouble(t[5]);
							//genes[ii[0]].get(ii[1]).sign_hotspot_syn[i]=Double.parseDouble(t[6]);
							//genes[ii[0]].get(ii[1]).sign_complete_syn[i]=Double.parseDouble(t[7]);ddd
							genes[ii[0]].get(ii[1]).sign_combined_uniform[i]=Double.parseDouble(t[8]);
							//genes[ii[0]].get(ii[1]).sign_destruct[i]=Double.parseDouble(t[9]);
							//genes[ii[0]].get(ii[1]).sign_hotspot[i]=Double.parseDouble(t[10]);
							//genes[ii[0]].get(ii[1]).sign_complete[i]=Double.parseDouble(t[11]);
						}
						
					}
					input.close();
				}
			}
			
			
			*/
			
			
			
			
			for (int k=0;k<entities.length;k++){
				int count_all=0;
				int count_destruct=0;
				in=new FileInputStream(file_destructive+entities[k]+".txt");
				inn=new DataInputStream(in);
				input= new BufferedReader(new InputStreamReader(inn));
				while((s=input.readLine())!=null){
					String[] t=s.split("	");
					int[] index=table_gene.get(t[0]);
					if(index!=null){//[0]!=-1
						//if(index[0]<genes.length&&index[1]<genes[index[0]].size()){
							genes[index[0]].get(index[1]).count_all[k]=Integer.parseInt(t[1]);
							genes[index[0]].get(index[1]).count_destruct[k]=Integer.parseInt(t[2]);
							count_all+=Integer.parseInt(t[1]);
							count_destruct+=Integer.parseInt(t[2]);
					
						//}
					}
				}
				input.close();
				
				double frac_destruct=(double)(count_destruct)/(double)(count_all);
				for (int i=0;i<genes.length;i++){
					for (int j=0;j<genes[i].size();j++){
						BinomialDistribution dist=new BinomialDistribution(genes[i].get(j).count_all[k],frac_destruct);
						if(genes[i].get(j).count_all[k]<2||genes[i].get(j).count_destruct[k]==0){
							genes[i].get(j).sign_destruct[k]=1-dist.cumulativeProbability(genes[i].get(j).count_destruct[k]-1);
						}
						else{
							genes[i].get(j).sign_destruct[k]=1-dist.cumulativeProbability(genes[i].get(j).count_destruct[k]-1);
						}
					}
				}
			}
			
			
			
			//test for accumulation of destructive mutations (just a very straightforward binomial distribution does the job)	
			
			
			//combination of the 3 p-values using brown
			
			for (int k=0;k<entities.length;k++){
				ArrayList<Double> product_syn=new ArrayList<Double>();
				for (int i=0;i<genes.length;i++){
					for (int j=0;j<genes[i].size();j++){
						double prod=genes[i].get(j).sign_vector_syn[k]*genes[i].get(j).sign_hotspot_syn[k];
						if(prod!=0&&prod!=1&&!Double.isNaN(prod)){
							product_syn.add(prod);
						}
					}
				}
				
				
				double avg_syn=0;
				for (int i=0;i<product_syn.size();i++){
					avg_syn+=(-2)*Math.log(product_syn.get(i));
				}
				avg_syn/=(double)(product_syn.size());
				
				double var_syn=0;
				for (int i=0;i<product_syn.size();i++){
					var_syn+=((-2)*Math.log(product_syn.get(i))-avg_syn)*((-2)*Math.log(product_syn.get(i))-avg_syn);
				}
				var_syn/=(double)(product_syn.size());
			
				double c1_syn=var_syn/(2*avg_syn);
				double k1_syn=2*avg_syn*avg_syn/var_syn;
				
				ChiSquaredDistribution dist1_syn=new ChiSquaredDistribution(k1_syn);
				for (int i=0;i<genes.length;i++){
					for (int j=0;j<genes[i].size();j++){
						double prod=genes[i].get(j).sign_vector_syn[k]*genes[i].get(j).sign_hotspot_syn[k];
						if(prod==1){
							genes[i].get(j).sign_complete_syn[k]=1;
						}
						else if(prod==0){
							genes[i].get(j).sign_complete_syn[k]=0;
						}
						else if(Double.isNaN(prod)){
							genes[i].get(j).sign_complete_syn[k]=Double.NaN;
						}
						else{
							genes[i].get(j).sign_complete_syn[k]=1-dist1_syn.cumulativeProbability((-2)*Math.log(prod)/c1_syn);
						}
					}
				}
				
			}
			
			
			
			//determine the right no. dimensions and scaling factor for the Brown method for the syn mutations
			
			
			//Combine p-values using brown for the syn mutations
			
			for (int k=0;k<entities.length;k++){
				ArrayList<Double> product=new ArrayList<Double>();
				for (int i=0;i<genes.length;i++){
					for (int j=0;j<genes[i].size();j++){
						double prod=genes[i].get(j).sign_combined[k]*Math.min(genes[i].get(j).sign_destruct[k],genes[i].get(j).sign_hotspot[k]);
						if(prod!=0&&prod!=1&&!Double.isNaN(prod)){
							product.add(prod);
						}
					}
				}
				
				double avg=0;
				for (int i=0;i<product.size();i++){
					avg+=(-2)*Math.log(product.get(i));
				}
				avg/=(double)(product.size());
				
				double var=0;
				for (int i=0;i<product.size();i++){
					var+=((-2)*Math.log(product.get(i))-avg)*((-2)*Math.log(product.get(i))-avg);
				}
				var/=(double)(product.size());
				
				
				//determine the right no. dimensions and scaling factor for the Brown method for the nonsyn mutations
				double c1=var/(2*avg);
				double k1=2*avg*avg/var;
				//Combine p-values using brown for the nonsyn mutations
				
				ChiSquaredDistribution dist1=new ChiSquaredDistribution(k1);
				for (int i=0;i<genes.length;i++){
					for (int j=0;j<genes[i].size();j++){
						double prod=genes[i].get(j).sign_combined[k]*Math.min(genes[i].get(j).sign_destruct[k],genes[i].get(j).sign_hotspot[k]);
						if(prod==1){
							genes[i].get(j).sign_complete[k]=1;
						}
						else if(prod==0){
							genes[i].get(j).sign_complete[k]=0;
						}
						else if(Double.isNaN(prod)){
							genes[i].get(j).sign_complete[k]=Double.NaN;
						}
						else{
							genes[i].get(j).sign_complete[k]=1-dist1.cumulativeProbability((-2)*Math.log(prod)/c1);
						}
					}
				}
			}
			
			for (int k=0;k<entities.length;k++){
				if(compute_uniform[k]){
					ArrayList<Double> product=new ArrayList<Double>();
					for (int i=0;i<genes.length;i++){
						for (int j=0;j<genes[i].size();j++){
							double prod=genes[i].get(j).sign_combined_uniform[k]*Math.min(genes[i].get(j).sign_destruct[k],genes[i].get(j).sign_hotspot[k]);
							if(prod!=0&&prod!=1&&!Double.isNaN(prod)){
								product.add(prod);
							}
						}
					}
					
					double avg=0;
					for (int i=0;i<product.size();i++){
						avg+=(-2)*Math.log(product.get(i));
					}
					avg/=(double)(product.size());
					
					double var=0;
					for (int i=0;i<product.size();i++){
						var+=((-2)*Math.log(product.get(i))-avg)*((-2)*Math.log(product.get(i))-avg);
					}
					var/=(double)(product.size());
					
					
					//determine the right no. dimensions and scaling factor for the Brown method for the nonsyn mutations
					double c1=var/(2*avg);
					double k1=2*avg*avg/var;
					//Combine p-values using brown for the nonsyn mutations
					
					ChiSquaredDistribution dist1=new ChiSquaredDistribution(k1);
					for (int i=0;i<genes.length;i++){
						for (int j=0;j<genes[i].size();j++){
							double prod=genes[i].get(j).sign_combined_uniform[k]*Math.min(genes[i].get(j).sign_destruct[k],genes[i].get(j).sign_hotspot[k]);
							if(prod==1){
								genes[i].get(j).sign_complete_uniform[k]=1;
							}
							else if(prod==0){
								genes[i].get(j).sign_complete_uniform[k]=0;
							}
							else if(Double.isNaN(prod)){
								genes[i].get(j).sign_complete_uniform[k]=Double.NaN;
							}
							else{
								genes[i].get(j).sign_complete_uniform[k]=1-dist1.cumulativeProbability((-2)*Math.log(prod)/c1);
							}
						}
					}
				}
				
			}
			
			Comparator<GeneSmall> comp_gene_syn=(GeneSmall g1, GeneSmall g2)->{
				return new Double(g1.sign_complete_syn).compareTo(new Double(g2.sign_complete_syn));
			};
			Comparator<GeneSmall> comp_gene=(GeneSmall g1, GeneSmall g2)->{
				return new Double(g1.sign_complete).compareTo(new Double(g2.sign_complete));
			};
			
			for (int k=0;k<entities.length;k++){
				int NN=0;
				for (int i=0;i<genes.length;i++){
					for (int j=0;j<genes[i].size();j++){
						if(genes[i].get(j).count_syn[k]+genes[i].get(j).count[k]+genes[i].get(j).count_all[k]+genes[i].get(j).count_destruct[k]>0){
							NN++;
						}
					}
				}
				
				ArrayList<GeneSmall> genes_all=new ArrayList<GeneSmall>();
				for (int i=0;i<genes.length;i++){
					for (int j=0;j<genes[i].size();j++){
						genes_all.add(genes[i].get(j).small(k));
					}
				}
				
				Collections.sort(genes_all,comp_gene_syn);
				for (int i=0;i<genes_all.size();i++){
					genes_all.get(i).fdr_syn=Math.min(1, (double)(genes_all.size())/(double)(i+1)*genes_all.get(i).sign_complete_syn);
				}
				
				//Compute fdr
				Collections.sort(genes_all,comp_gene);
				for (int i=0;i<genes_all.size();i++){
					if(i+1<=NN){
						genes_all.get(i).fdr=Math.min(1, (double)(NN)/(double)(i+1)*genes_all.get(i).sign_complete);
					}
					else{
						//genes_all.get(i).fdr=1;
						break;
					}
					
				}
			
				FileWriter out=new FileWriter(file_out+entities[k]+".txt");
				BufferedWriter output= new BufferedWriter(out);
				//SignVector	SignCount	
				//output.write("Name	TargetSize	TargetSizeSyn	Count	CountSyn	SignVectorSyn	SignCumSyn	SignCompleteSyn	SignSeq	SignDm	SignCum	SignCombined	FDRSyn	FDR");
				output.write("Name	TargetSize	TargetSizeSyn	Count	CountSyn	SignificanceSyn	FDRSyn	Significance	FDR");
				output.newLine();
				
				for (int i=0;i<genes_all.size();i++){//genes_all.get(i).sign_vector+"	"+genes_all.get(i).sign_cbase+"	"+
					//output.write(genes_all.get(i).name+"	"+genes_all.get(i).cov+"	"+genes_all.get(i).cov_syn+"	"+genes_all.get(i).count+"	"+genes_all.get(i).count_syn+"	"+genes_all.get(i).sign_vector_syn+"	"+genes_all.get(i).sign_hotspot_syn+"	"+genes_all.get(i).sign_complete_syn+"	"+genes_all.get(i).sign_combined+"	"+genes_all.get(i).sign_destruct+"	"+genes_all.get(i).sign_hotspot+"	"+genes_all.get(i).sign_complete+"	"+genes_all.get(i).fdr_syn+"	"+genes_all.get(i).fdr);
					output.write(genes_all.get(i).name+"	"+genes_all.get(i).cov+"	"+genes_all.get(i).cov_syn+"	"+genes_all.get(i).count+"	"+genes_all.get(i).count_syn+"	"+genes_all.get(i).sign_complete_syn+"	"+genes_all.get(i).fdr_syn+"	"+genes_all.get(i).sign_complete+"	"+genes_all.get(i).fdr);
					output.newLine();
				}
				output.close();
				
				
				if(compute_uniform[k]){

					ArrayList<GeneSmall> genes_all_uniform=new ArrayList<GeneSmall>();
					for (int i=0;i<genes.length;i++){
						for (int j=0;j<genes[i].size();j++){
							genes_all_uniform.add(genes[i].get(j).small_uniform(k));
						}
					}
				
					Collections.sort(genes_all_uniform,comp_gene_syn);
					for (int i=0;i<genes_all_uniform.size();i++){
						genes_all_uniform.get(i).fdr_syn=Math.min(1, (double)(genes_all_uniform.size())/(double)(i+1)*genes_all_uniform.get(i).sign_complete_syn);
					}
					
					//Compute fdr
					Collections.sort(genes_all_uniform,comp_gene);
					for (int i=0;i<genes_all_uniform.size();i++){
						if(i+1<=NN){
							genes_all_uniform.get(i).fdr=Math.min(1, (double)(NN)/(double)(i+1)*genes_all_uniform.get(i).sign_complete);
						}
						else{
							//genes_all_uniform.get(i).fdr=1;
							break;
						}
						
					}
			
					out=new FileWriter(file_out_uniform+entities[k]+".txt");
					output= new BufferedWriter(out);
					//SignVector	SignCount	
					//output.write("Name	TargetSize	TargetSizeSyn	Count	CountSyn	SignVectorSyn	SignCumSyn	SignSeqSyn	SignSeq	SignDm	SignCum	SignCombined	FDRSyn	FDR");
					output.write("Name	TargetSize	TargetSizeSyn	Count	CountSyn	SignificanceSyn	FDRSyn	Significance	FDR");
					output.newLine();
					
					for (int i=0;i<genes_all_uniform.size();i++){//genes_all.get(i).sign_vector+"	"+genes_all.get(i).sign_cbase+"	"+
						//output.write(genes_all_uniform.get(i).name+"	"+genes_all_uniform.get(i).cov+"	"+genes_all_uniform.get(i).cov_syn+"	"+genes_all_uniform.get(i).count+"	"+genes_all_uniform.get(i).count_syn+"	"+genes_all_uniform.get(i).sign_vector_syn+"	"+genes_all_uniform.get(i).sign_hotspot_syn+"	"+genes_all_uniform.get(i).sign_complete_syn+"	"+genes_all_uniform.get(i).sign_combined+"	"+genes_all_uniform.get(i).sign_destruct+"	"+genes_all_uniform.get(i).sign_hotspot+"	"+genes_all_uniform.get(i).sign_complete+"	"+genes_all_uniform.get(i).fdr_syn+"	"+genes_all_uniform.get(i).fdr);
						output.write(genes_all_uniform.get(i).name+"	"+genes_all_uniform.get(i).cov+"	"+genes_all_uniform.get(i).cov_syn+"	"+genes_all_uniform.get(i).count+"	"+genes_all_uniform.get(i).count_syn+"	"+genes_all_uniform.get(i).sign_complete_syn+"	"+genes_all_uniform.get(i).fdr_syn+"	"+genes_all_uniform.get(i).sign_complete+"	"+genes_all_uniform.get(i).fdr);
						output.newLine();
					}
					output.close();
				}
			}
			
			
			//for FDR correction put all genes in 1 list and sort them by sign (both for syn sign and non syn sign)
//			ArrayList<Gene> genes_all=new ArrayList<Gene>();
//			for (int i=0;i<genes.length;i++){
//				genes_all.addAll(genes[i]);
//			}
			
			
			
			
			
			
			
			
			
			/*
			int NN=0;
			for (int i=0;i<genes_all.size();i++){
				if(genes_all.get(i).count_syn+genes_all.get(i).count+genes_all.get(i).count_all+genes_all.get(i).count_destruct>0){
					NN++;
				}
			}*/
			
			
			
			//output of FDR values and 
			
			
			
				
				
			/*
			run(index("1",chr));
			System.out.println("END");
			int ii=index("Skin",entities);
//			
			for (int i=0;i<pos2.size();i++){
				System.out.println(pos2.get(i)+"	"+nucl2.get(i)+"	"+coverage2.get(i)+"	"+label2.get(i)[0]+";"+label2.get(i)[1]+";"+label2.get(i)[2]+"	"+lambda2.get(i)[ii][0]+";"+lambda2.get(i)[ii][1]+";"+lambda2.get(i)[ii][2]+"	"+count2.get(i)[ii][0]+";"+count2.get(i)[ii][1]+";"+count2.get(i)[ii][2]);
			}
			*/
			
			
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
			
			System.out.println(e);
		}
	}
	
	public static boolean isDouble(String x){
		try{
			Double.parseDouble(x);
			return true;
		}
		catch(Exception e){
			return false;
		}
	}
	
	public static double log_gamma_int(int x){
		if(x<gamma_precompute.length){
			return gamma_precompute[x];
		}
		else{
			return Gamma.logGamma(x);
		}
	}
	
	//each thread computes the clusterwise local mutation rate for 1 chromosome
	//(reading the reference sequence is usually the rate limiting step)
	//private static class Subthread extends Thread{
	
	static ArrayList<String> nucl2=new ArrayList<String>();
	static ArrayList<Integer> pos2=new ArrayList<Integer>();
	static ArrayList<double[][]> lambda2=new ArrayList<double[][]>();
	//static ArrayList<double[][]> lambda22=new ArrayList<double[][]>();
	static ArrayList<Double> coverage2=new ArrayList<Double>();
	static ArrayList<int[]> label2=new ArrayList<int[]>();
	static ArrayList<int[][]> count2=new ArrayList<int[][]>();
	static ArrayList<String> amino_acid2=new ArrayList<String>();
	
		
		//BufferedWriter output=null;
		public static void run(int c){
			
//			for (int nn=0;nn<genes[c].size();nn++){
//				System.out.println(nn+"	"+genes[c].get(nn).start+"	"+genes[c].get(nn).end);
//			}
//			System.exit(0);
			//int center=50;
			/*
			for (int i=0;i<lambda_type.size();i++){
				for (int j=0;j<lambda_type.get(i).length;j++){
					for (int k=0;k<lambda_type.get(i)[j].length;k++){
						lambda_type.get(i)[j][k]=(double)((int)(1000*lambda_type.get(i)[j][k]))/1000.0;
						
					}
				}
			}
			
			for (int i=0;i<lambda_context.size();i++){
				for (int j=0;j<lambda_context.get(i).length;j++){
					for (int k=0;k<lambda_context.get(i)[j].length;k++){
						for (int l=0;l<lambda_context.get(i)[j][k].length;l++){
							lambda_context.get(i)[j][k][l]=(double)((int)(1000*lambda_context.get(i)[j][k][l]))/1000.0;
						}
					}
				}
			}*/
			/*
			int n1=0;
			int n2=0;
			for (int i=0;i<lambda_context.size();i++){
				for (int j=0;j<lambda_context.get(i).length;j++){
					for (int k=0;k<lambda_context.get(i)[j].length;k++){
						for (int l=0;l<lambda_context.get(i)[j][k].length;l++){
							double x=Math.max(1.0/lambda_context.get(i)[j][k][l], lambda_context.get(i)[j][k][l]);
							if(x>=1.1){
								n1++;
							}
							else{
								n2++;
							}
						}
					}
				}
			}
			System.out.println(n1+"	"+n2);
			*/
			
			
			pos=new ArrayList<Integer>();
			nucl=new ArrayList<String>();
			nucl_index=new ArrayList<Integer>();
			coverage=new ArrayList<Double>();
			label=new ArrayList<int[]>();
			count=new ArrayList<int[][]>();
			amino_acid=new ArrayList<String>();
			
			nucl2=new ArrayList<String>();
			pos2=new ArrayList<Integer>();
			lambda2=new ArrayList<double[][]>();
			//lambda22=new ArrayList<double[][]>();
			coverage2=new ArrayList<Double>();
			label2=new ArrayList<int[]>();
			count2=new ArrayList<int[][]>();
			amino_acid2=new ArrayList<String>();
			
			System.gc();
			System.runFinalization();
			
			try{
				String s="";
				for (int i=0;i<50;i++){
					pos.add(-1);
					nucl.add("N");
					nucl_index.add(-1);
					coverage.add(0.0);
					label.add(new int[0]);
					count.add(new int[0][0]);
					amino_acid.add("-");
				}
				
				//walk through the refernce sequence and safe it into a queue
				//FileWriter out=new FileWriter(file_out+chr[c]+".txt");
				//output= new BufferedWriter(out);
				
				int n=0;
				{
					FileInputStream in=new FileInputStream(file_annotation+chr[c]+".txt");
					DataInputStream inn=new DataInputStream(in);
					BufferedReader input= new BufferedReader(new InputStreamReader(inn));
					while(input.readLine()!=null){
						n++;
					}
					input.close();
				}
				
				FileInputStream in=new FileInputStream(file_annotation+chr[c]+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				
				FileInputStream in2=new FileInputStream(file_align+chr[c]+".txt");
				DataInputStream inn2=new DataInputStream(in2);
				BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
				
				int nnn=0;
				for (int nn=0;nn<genes[c].size();nn++){
					
					int ii=0;
					if(pos2.size()>0){
						while(ii<pos2.size()&&pos2.get(ii)<genes[c].get(nn).start){
							ii++;
						}
					}
					
					//remove the squence of the previous gene out of the queue
					if(ii>0){
						for (int i=ii-1;i>=0;i--){
							nucl2.remove(i);
							lambda2.remove(i);
							//lambda22.remove(i);
							pos2.remove(i);
							coverage2.remove(i);
							count2.remove(i);
							label2.remove(i);
							amino_acid2.remove(i);
						}
					}
					
					
					while((s=input.readLine())!=null){
						nnn++;
						//if(nnn%10000==0){
						//	System.out.println(nnn+"/"+n+"	"+System.currentTimeMillis());
						//}
						
						String[] t=s.split("	");
						String[] t2=input2.readLine().split("	");
						
						int[] label_local=new int[3];
						if(t.length>3){
							if(t[3].equals(t[5])){
								label_local[0]=0;//"syn";
							}
							else{
								label_local[0]=1;//"ns";
							}
						}
						else{
							label_local[0]=2;//"nc";
						}
						
						if(t.length>3){
							if(t[3].equals(t[6])){
								label_local[1]=0;//"syn";
							}
							else{
								label_local[1]=1;//"ns";
							}
						}
						else{
							label_local[1]=2;//"nc";
						}
						
						if(t.length>3){
							if(t[3].equals(t[7])){
								label_local[2]=0;//"syn";
							}
							else{
								label_local[2]=1;//"ns";
							}
						}
						else{
							label_local[2]=2;//"nc";
						}
						
						int[][] count_local=new int[entities.length+1][3];
						for (int k=0;k<3;k++){
							if(t2.length>k&&!t2[k].equals("")){
								String[] tt=t2[k].split(";");
								for (int l=0;l<tt.length;l++){
									count_local[table_entity.get(Integer.parseInt(tt[l]))][k]++;
								}
								count_local[count_local.length-1][k]=tt.length;
							}
							
						}
	
						
						pos.add(Integer.parseInt(t[0]));
						nucl.add(t[1]);
						nucl_index.add(nucl_index(t[1]));
						coverage.add(Double.parseDouble(t[2]));
						label.add(label_local);
						count.add(count_local);
						if(t.length>3){
							amino_acid.add(t[3]+t[4]);
						}
						else{
							amino_acid.add("-");
						}
						
						
						//when the queue is large enough, compute the mutation rate in the cnenter
						//(update function) and delete the first element
						if(pos.size()>100){
							//update_10_10();
							update_10_10();
							/*
							lambda2.add(new double[weights.length][3]);
							pos2.add(pos.get(50));
							nucl2.add(nucl.get(50));
							coverage2.add(coverage.get(50));
							label2.add(label.get(50));
							count2.add(count.get(50));
							amino_acid2.add(amino_acid.get(50));
							*/
							
							
							/*
							System.out.print(pos2.get(pos2.size()-1)+"	"+nucl2.get(nucl2.size()-1)+"	"+coverage2.get(coverage2.size()-1));
							for (int j=0;j<3;j++){
								System.out.print("	"+lambda2.get(lambda2.size()-1)[iii][j]);
							}
							System.out.println();
							*/
							
							pos.remove(0);
							nucl.remove(0);
							nucl_index.remove(0);
							coverage.remove(0);
							label.remove(0);
							count.remove(0);
							amino_acid.remove(0);
							//System.out.println(pos2.get(pos2.size()-1)+"	"+genes[c].get(nn).end);
							if(pos2.get(pos2.size()-1)>=genes[c].get(nn).end){
								break;
							}
							
							
						}
					}
					if(s==null){
						for (int i=0;i<50;i++){
							pos.add(-1);
							nucl.add("N");
							nucl_index.add(-1);
							coverage.add(0.0);
							label.add(new int[0]);
							count.add(new int[0][0]);
							amino_acid.add("-");
							//update_10_10();
							update_10_10();
							pos.remove(0);
							nucl.remove(0);
							nucl_index.remove(0);
							coverage.remove(0);
							label.remove(0);
							count.remove(0);
							amino_acid.remove(0);
						}
					}
					
					/*
					if(table_failed_genes.get(genes[c].get(nn).name)==null){
						continue;
					}*/
					
					//ArrayList<String> aa_gene=new ArrayList<String>();
					//ArrayList<String> aa_gene_syn=new ArrayList<String>();
					ArrayList<int[]> index_gene=new ArrayList<int[]>();
					ArrayList<int[]> index_gene_syn=new ArrayList<int[]>();
					ArrayList<int[]> index_gene2=new ArrayList<int[]>();
					ArrayList<int[]> index_gene_syn2=new ArrayList<int[]>();
					
					for (int i=0;i<pos2.size();i++){
						if(!genes[c].get(nn).contains(pos2.get(i))){
							continue;
						}
						if(coverage2.get(i)>0.5){
							for (int j=0;j<3;j++){
								if(label2.get(i)[j]==1){
									index_gene.add(new int[]{i,j});
									//aa_gene.add(amino_acid2.get(i));
								}
								else if(label2.get(i)[j]==0){
									index_gene_syn.add(new int[]{i,j});
									//aa_gene_syn.add(amino_acid2.get(i));
								}
							}
						}
						
						if(coverage2.get(i)>=0.05){//0.3
							for (int j=0;j<3;j++){
								if(label2.get(i)[j]!=0){
									index_gene2.add(new int[]{i,j});
									//lambda_gene.add(lambda2.get(i)[k][j]);
									//count_gene.add(count2.get(i)[k][j]);
								}
								else{
									index_gene_syn2.add(new int[]{i,j});
									//lambda_gene_syn.add(lambda2.get(i)[k][j]);
									//count_gene_syn.add(count2.get(i)[k][j]);
								}
							}
							
							
						}
					}
					
					boolean is_or=false;
					if(genes[c].get(nn).name.length()>2&&genes[c].get(nn).name.substring(0, 2).equals("OR")&&contains(genes[c].get(nn).name.substring(2,3),numbers)){//&&genes[i].get(j).ell_s>0.0
						is_or=true;
					}
						
					for (int k=0;k<entities.length;k++){
//						ArrayList<double[]> lambda22=new ArrayList<double[]>();
//						for (int i=0;i<pos2.size();i++){
//							if(lambda2.get(i).length==1){
//								int iii=(int)(lambda2.get(i)[0][0]);
//								if(nucl2.get(i).equals("C")||nucl2.get(i).equals("G")){
//									lambda22.add(new double[]{lambda_context_product6_weight[k][tt2[0]][iii],lambda_context_product6_weight[k][tt2[1]][iii],lambda_context_product6_weight[k][tt2[2]][iii]});
//								}
//								else {
//									lambda22.add(new double[]{lambda_context_product6_weight[k][tt1[0]][iii],lambda_context_product6_weight[k][tt1[1]][iii],lambda_context_product6_weight[k][tt1[2]][iii]});
//								}
//							}
//							else{
//								lambda22.add(lambda2.get(i)[k]);
//							}
//						}
						
						
						//if(!entities[k].equals("Skin")){
						//	continue;
						//}
						//System.out.println(k+"	"+entities[k]);
						//todo: does this include pan-cancer? if so, remove...
						//if(!sign_entities[k]&&!sign_uniform_entities[k]){//TODO: select whether or not you want to compute sign or sign uniform
						//	continue;
						//}
						
						//ArrayList<Double> lambda_gene=new ArrayList<Double>();
						//ArrayList<Integer> count_gene=new ArrayList<Integer>();
						//ArrayList<Double> coverage_gene=new ArrayList<Double>();
						//ArrayList<Double> lambda_gene_syn=new ArrayList<Double>();
						//ArrayList<Integer> count_gene_syn=new ArrayList<Integer>();
						//ArrayList<Double> coverage_gene_syn=new ArrayList<Double>();
						//ArrayList<double[]> pos_gene=new ArrayList<double[]>();
						//ArrayList<double[]> pos_gene_syn=new ArrayList<double[]>();
						//ArrayList<String> aa_gene=new ArrayList<String>();
						//ArrayList<String> aa_gene_syn=new ArrayList<String>();
						
//						for (int i=0;i<pos2.size();i++){
//							if(!genes[c].get(nn).contains(pos2.get(i))){
//								continue;
//							}
//							if(lambda2.get(i).length==1){
//								int iii=(int)(lambda2.get(i)[0][0]);
//								//System.out.println(lambda2.get(i)[0][0]+"	"+iii);
//								
//								
////								if(entities[k].equals("Skin")){
////									if(nucl2.get(i).equals("C")||nucl2.get(i).equals("G")){
////										System.out.print(pos2.get(i)+"	"+nucl2.get(i));
////										for (int j=0;j<3;j++){
////											System.out.print("	"+lambda_context_product6_weight[k][tt2[j]][iii]);
////										}
////										for (int j=0;j<3;j++){
////											System.out.print("	"+lambda22.get(i)[k][j]);
////										}
////										System.out.println();
////										
////									}
////									else if(nucl2.get(i).equals("A")||nucl2.get(i).equals("T")){
////										System.out.print(pos2.get(i)+"	"+nucl2.get(i));
////										for (int j=0;j<3;j++){
////											System.out.print("	"+lambda_context_product6_weight[k][tt1[j]][iii]);
////										}
////										for (int j=0;j<3;j++){
////											System.out.print("	"+lambda22.get(i)[k][j]);
////										}
////										System.out.println();
////									}
////								}
//								
//								
//								if(nucl2.get(i).equals("C")||nucl2.get(i).equals("G")){
//									
//									for (int j=0;j<3;j++){
//										if(label2.get(i)[j]==1&&coverage2.get(i)>0.5){
//											pos_gene.add(new double[]{count2.get(i)[k][j],lambda_context_product6_weight[k][tt2[j]][iii]});
//											//aa_gene.add(amino_acid2.get(i));
//										}
//										else if(label2.get(i)[j]==0&&coverage2.get(i)>0.5){
//											pos_gene_syn.add(new double[]{count2.get(i)[k][j],lambda_context_product6_weight[k][tt2[j]][iii]});
//											//aa_gene_syn.add(amino_acid2.get(i));
//										}
//									}
//									
//									if(coverage2.get(i)>=0.05){//0.3
//										for (int j=0;j<3;j++){
//											if(label2.get(i)[j]!=0){
//												lambda_gene.add(lambda_context_product6_weight[k][tt2[j]][iii]);//lambda2.get(i)[k][j]
//												count_gene.add(count2.get(i)[k][j]);
//											}
//											else{
//												lambda_gene_syn.add(lambda_context_product6_weight[k][tt2[j]][iii]);//lambda2.get(i)[k][j]
//												count_gene_syn.add(count2.get(i)[k][j]);
//											}
//										}
//									}
//									
//								}
//								else if(nucl2.get(i).equals("A")||nucl2.get(i).equals("T")){
//									
//									for (int j=0;j<3;j++){
//										if(label2.get(i)[j]==1&&coverage2.get(i)>0.5){
//											pos_gene.add(new double[]{count2.get(i)[k][j],lambda_context_product6_weight[k][tt1[j]][iii]});
//											//aa_gene.add(amino_acid2.get(i));
//										}
//										else if(label2.get(i)[j]==0&&coverage2.get(i)>0.5){
//											pos_gene_syn.add(new double[]{count2.get(i)[k][j],lambda_context_product6_weight[k][tt1[j]][iii]});
//											//aa_gene_syn.add(amino_acid2.get(i));
//										}
//									}
//									
//									if(coverage2.get(i)>=0.05){//0.3
//										for (int j=0;j<3;j++){
//											if(label2.get(i)[j]!=0){
//												lambda_gene.add(lambda_context_product6_weight[k][tt1[j]][iii]);//lambda2.get(i)[k][j]
//												count_gene.add(count2.get(i)[k][j]);
//											}
//											else{
//												lambda_gene_syn.add(lambda_context_product6_weight[k][tt1[j]][iii]);//lambda2.get(i)[k][j]
//												count_gene_syn.add(count2.get(i)[k][j]);
//											}
//										}
//									}
//								}
//								
//								
//							}
//							else{
//								
//								for (int j=0;j<3;j++){
//									if(label2.get(i)[j]==1&&coverage2.get(i)>0.5){
//										pos_gene.add(new double[]{count2.get(i)[k][j],lambda2.get(i)[k][j]});
//										//aa_gene.add(amino_acid2.get(i));
//									}
//									else if(label2.get(i)[j]==0&&coverage2.get(i)>0.5){
//										pos_gene_syn.add(new double[]{count2.get(i)[k][j],lambda2.get(i)[k][j]});
//										//aa_gene_syn.add(amino_acid2.get(i));
//									}
//								}
//								
//								if(coverage2.get(i)>=0.05){//0.3
//										
//									if(label2.get(i)[0]!=0){
//										lambda_gene.add(lambda2.get(i)[k][0]);
//										count_gene.add(count2.get(i)[k][0]);
//									}
//									else{
//										lambda_gene_syn.add(lambda2.get(i)[k][0]);
//										count_gene_syn.add(count2.get(i)[k][0]);
//									}
//									if(label2.get(i)[1]!=0){
//										lambda_gene.add(lambda2.get(i)[k][1]);
//										count_gene.add(count2.get(i)[k][1]);
//									}
//									else{
//										lambda_gene_syn.add(lambda2.get(i)[k][1]);
//										count_gene_syn.add(count2.get(i)[k][1]);
//									}
//									if(label2.get(i)[2]!=0){
//										lambda_gene.add(lambda2.get(i)[k][2]);
//										count_gene.add(count2.get(i)[k][2]);
//									}
//									else{
//										lambda_gene_syn.add(lambda2.get(i)[k][2]);
//										count_gene_syn.add(count2.get(i)[k][2]);
//									}
//								}
//							}
//							
//							
//							
//							

						
						ArrayList<double[]> lambda_count=new ArrayList<double[]>();
						for (int jj=0;jj<index_gene2.size();jj++){
							double x=0;
							if(lambda2.get(index_gene2.get(jj)[0]).length==1){
								int iii=(int)(lambda2.get(index_gene2.get(jj)[0])[0][0]);
								if(nucl2.get(index_gene2.get(jj)[0]).equals("C")||nucl2.get(index_gene2.get(jj)[0]).equals("G")){
									x=lambda_context_product6_weight[k][tt2[index_gene2.get(jj)[1]]][iii];
								}
								else {
									x=lambda_context_product6_weight[k][tt1[index_gene2.get(jj)[1]]][iii];
								}
							}
							else{
								x=lambda2.get(index_gene2.get(jj)[0])[k][index_gene2.get(jj)[1]];
							}
							
							lambda_count.add(new double[]{x,count2.get(index_gene2.get(jj)[0])[k][index_gene2.get(jj)[1]]});//lambda22.get(index_gene2.get(i)[0])[index_gene2.get(i)[1]]
						}
						Collections.sort(lambda_count,comppp);
						
						ArrayList<double[]> lambda_count_syn=new ArrayList<double[]>();
						for (int jj=0;jj<index_gene_syn2.size();jj++){
							double x=0;
							if(lambda2.get(index_gene_syn2.get(jj)[0]).length==1){
								int iii=(int)(lambda2.get(index_gene_syn2.get(jj)[0])[0][0]);
								if(nucl2.get(index_gene_syn2.get(jj)[0]).equals("C")||nucl2.get(index_gene_syn2.get(jj)[0]).equals("G")){
									x=lambda_context_product6_weight[k][tt2[index_gene_syn2.get(jj)[1]]][iii];
								}
								else {
									x=lambda_context_product6_weight[k][tt1[index_gene_syn2.get(jj)[1]]][iii];
								}
							}
							else{
								x=lambda2.get(index_gene_syn2.get(jj)[0])[k][index_gene_syn2.get(jj)[1]];
							}
							
							lambda_count_syn.add(new double[]{x,count2.get(index_gene_syn2.get(jj)[0])[k][index_gene_syn2.get(jj)[1]]});//lambda22.get(index_gene_syn2.get(i)[0])[index_gene_syn2.get(i)[1]]
						}
						Collections.sort(lambda_count_syn,comppp);
						
						
//						ArrayList<double[]> x=sort(lambda_gene,count_gene);
//						
//						lambda_gene=new ArrayList<Double>();
//						count_gene=new ArrayList<Integer>();
//						for (int i=0;i<x.size();i++){
//							lambda_gene.add(x.get(i)[0]);
//							count_gene.add((int)(x.get(i)[1]));
//						}
//						x=sort(lambda_gene_syn,count_gene_syn);
//						
//						lambda_gene_syn=new ArrayList<Double>();
//						count_gene_syn=new ArrayList<Integer>();
//						for (int i=0;i<x.size();i++){
//							lambda_gene_syn.add(x.get(i)[0]);
//							count_gene_syn.add((int)(x.get(i)[1]));
//						}
						
						int err=3;//TODO: 3 ??
						
						//System.out.println(genes[c].get(nn).name+"	"+genes[c].get(nn).cov+"	"+genes[c].get(nn).cov_syn);
						
						if(genes[c].get(nn).cov*genes[c].get(nn).cov_syn==0){
							//genes[c].get(nn).sign_vector[k]=1;
							genes[c].get(nn).sign_hotspot[k]=1;//signREVISED(lambda2,count2,genes[c].get(nn).prob_nonsyn,err);//err,
							genes[c].get(nn).sign_combined[k]=1;//err,
							//genes[c].get(nn).sign_cbase[k]=1;
							genes[c].get(nn).sign_vector_syn[k]=1;
							genes[c].get(nn).sign_hotspot_syn[k]=1;
						
						}
						
						else{
							
							//if(entities[k].equals("Skin")){
								
								//System.out.println("Sign "+genes[c].get(nn).name+"	"+entities[k]);
								
								//TODO: reevaluate which of these sign values are needed later
								
								
								//System.out.println("size	"+prob_nonsyn.size());
								
								
								genes[c].get(nn).sign_hotspot[k]=hotspot_sign(summarize(count2,lambda2,nucl2,amino_acid2, index_gene,k));//signREVISED(lambda2,count2,genes[c].get(nn).prob_nonsyn,err);//err,
								
								
								//genes[c].get(nn).sign_hotspot[k]=hotspot_sign(summarize(pos_gene,aa_gene));//signREVISED(lambda2,count2,genes[c].get(nn).prob_nonsyn,err);//err,
								if(!is_or){
									//if(entities[k].equals("Skin")){
										
										/*
										double sobs = genes[c].get(nn).sobs[k];//["sobs"];
										double xobs = genes[c].get(nn).xobs[k];//["xobs"];
										double sexp = genes[c].get(nn).ell_s[k];//["ell_s"];
										double xexp = genes[c].get(nn).ell_x[k];
										System.out.println(genes[c].get(nn).name+"	"+sobs+"	"+xobs+"	"+sexp+"	"+xexp);
										*/
										
										//for (int l=0;l<prob_nonsyn.size();l++){
										//	System.out.println(l+"	"+prob_nonsyn.get(l));
										//}
										ArrayList<Double> prob_nonsyn=prob_nonsyn(genes[c].get(nn),params[k],mod_C[k], k);
									
										if(Double.isNaN(prob_nonsyn.get(prob_nonsyn.size()-1))||Double.isInfinite(prob_nonsyn.get(prob_nonsyn.size()-1))){
											genes[c].get(nn).sign_combined[k]=1.0;
										}
										else{
											genes[c].get(nn).sign_combined[k]=signREVISED(genes[c].get(nn).count[k],lambda_count,prob_nonsyn,err);//signREVISED(genes[c].get(nn).count[k],lambda_gene,count_gene,prob_nonsyn,err);//genes[c].get(nn). //err,
										}
										
										if(compute_uniform[k]){
											
											ArrayList<Double> prob_nonsyn_uniform=new ArrayList<Double>();
											double prob=1;
											int nnnn=0;
											while(prob>=Math.pow(10, -6)||nnnn<=genes[c].get(nn).count[k]){
												prob=probability_binom(non_syn_count[k],nnnn,genes[c].get(nn).cov/non_syn_coverage);
												//if(entities[k].equals("Thyroid")){
												//	System.out.println(non_syn_count[k]+"	"+nnnn+"	"+genes[c].get(nn).cov/non_syn_coverage);
												//}
												prob_nonsyn_uniform.add(prob);
												nnnn++;
											}
											
//											if(entities[k].equals("Thyroid")){
//												System.out.println(genes[c].get(nn).name);
//												for (int i=0;i<prob_nonsyn_uniform.size();i++){
//													System.out.println(i+"	"+prob_nonsyn_uniform.get(i));
//												}
//											}
											
											if(Double.isNaN(prob_nonsyn_uniform.get(prob_nonsyn_uniform.size()-1))||Double.isInfinite(prob_nonsyn_uniform.get(prob_nonsyn_uniform.size()-1))){
												genes[c].get(nn).sign_combined_uniform[k]=1.0;
											}
											else{
												genes[c].get(nn).sign_combined_uniform[k]=signREVISED(genes[c].get(nn).count[k],lambda_count,prob_nonsyn_uniform,err);//signREVISED(genes[c].get(nn).count[k],lambda_gene,count_gene,prob_nonsyn,err);//genes[c].get(nn). //err,
												/*if(entities[k].equals("Thyroid")){
													System.out.println(genes[c].get(nn).name+"	"+genes[c].get(nn).sign_combined_uniform[k]);
												}*/
											}
											
										}
										
										
									//}
								}
								else{
									genes[c].get(nn).sign_combined[k]=1.0;
								}
								genes[c].get(nn).sign_vector_syn[k]=sign(lambda_count_syn,err);//sign(lambda_gene_syn,count_gene_syn,err);
								//genes[c].get(nn).sign_hotspot_syn[k]=hotspot_sign(summarize(pos_gene_syn,aa_gene_syn));
								genes[c].get(nn).sign_hotspot_syn[k]=1.0;//hotspot_sign(summarize(count2,lambda2,nucl2,amino_acid2, index_gene_syn,k));//signREVISED(lambda2,count2,genes[c].get(nn).prob_nonsyn,err);//err,
							
								//if(entities[k].equals("Skin")){
								//		System.out.println(genes[c].get(nn).name+"	"+genes[c].get(nn).count[k]+"	"+genes[c].get(nn).sign_combined[k]);
									//genes[c].get(nn).sign_hotspot[k]+"	"++"	"+genes[c].get(nn).sign_hotspot_syn[k]+"	"+genes[c].get(nn).sign_vector_syn[k]
								//}
								
						}
							
					}
				
				

					
					
				}
				//compute the last mutation rates at the end of the chr until the queue empty 
				
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
		
			System.out.println("Done "+chr[c]);
		}
		
		public static double probability_binom(int n , int k, double f){
			return Math.exp(Gamma.logGamma(n+1)-Gamma.logGamma(k+1)-Gamma.logGamma(n-k+1)+Math.log(f)*k+Math.log(1-f)*(n-k));
		}
		
		public static boolean contains(String s, String[] t){
			for (int i=0;i<t.length;i++){
				if(t[i].equals(s)){
					return true;
				}
			}
			return false;
		}
		
		public static double[] sub(double[] x, int start, int end){
			double[] y=new double[end-start];
			for (int i=start;i<end;i++){
				y[i-start]=x[i];
			}
			return y;
		}
		
		public static ArrayList<Double> prob_nonsyn(Gene g , double[] p, int modC, int k){
			ArrayList<Double> probs=new ArrayList<Double>();
			
			double L=1.0;
			double sobs = g.sobs[k];//["sobs"];
			double xobs = g.xobs[k];//["xobs"];
			double sexp = g.ell_s[k];//["ell_s"];
			double xexp = g.ell_x[k];//["ell_x"];
			
			double ratx = xexp/sexp;
			//System.out.println(g.name+"	"+sobs+"	"+xobs+"	"+sexp+"	"+xexp+"	"+ratx);
			
			if (ratx<1e-10){
				return probs;
			}

			double last_p = 0.0;
			boolean flag_peak = false;
			for (int x=0;x<=10000000;x++){
				if (last_p>1.0 || Double.isNaN(last_p)){
					break;
				}
				double cur_p=0;
				if (last_p>1e-5){
						cur_p = pofx_given_s(x, sobs, L, ratx, p, modC);
						probs.add(cur_p);//new double[]{x, cur_p}
					
				}
				else{
					if (flag_peak && x>xobs){
						break;
					}
					else{
						cur_p = pofx_given_s(x, sobs, L, ratx,  p, modC);
						probs.add(cur_p);
						
					}
				}
				if (cur_p-last_p<0){
					flag_peak = true;
				}
				last_p = cur_p;
			}
			
			//if(Double.isNaN(probs.get(probs.size()-1))){
			//	probs.remove(probs.size()-1);
			//}
			
			return probs;
		}
		
		public static double sign(ArrayList<double[]> lambda_count ,int err){
			double sign_quick=sign_quick(lambda_count, err);
			if(sign_quick<0.001){//0.05
				return sign_long(lambda_count, err);
			}
			else{
				return sign_quick;
			}
			//return sign_quick;
		}
		
		static Comparator<double[]> comp0=(double[] a, double[] b)->{
			return -new Double(a[0]).compareTo(b[0]);
		};
		
		public static double sign_long(ArrayList<double[]> lambda_count, int err){
			
			
			int NN=0;
			for (int i=0;i<lambda_count.size();i++){
				NN+=lambda_count.get(i)[1];
			}
			if(lambda_count.size()==0||NN==0){
				return 1;
			}
			
			double T=0;
			for (int i=0;i<lambda_count.size();i++){
				if(lambda_count.get(i)[1]>0){
					T+=lambda_count.get(i)[1]*Math.log(lambda_count.get(i)[0]);
					//System.out.println(Math.log(lambda2.get(i)));
				}
				if(lambda_count.get(i)[1]>0){
					T-=log_gamma_int((int)(lambda_count.get(i)[1])+1);
				}
			}
			double offset=NN*Math.log((double)(1)/(double)(lambda_count.size()));
			T+=offset;
			//System.out.println(T);
			
			//System.out.println("T: "+Math.exp(T));
			
			//int NN=500;
			
			Collections.sort(lambda_count,comp0);//sort by 0 according to comp
			double sum_lambda=0;
			for (int i=0;i<lambda_count.size();i++){
				sum_lambda+=lambda_count.get(i)[0];
			}
			
			double sum_pow=0;
			{
				int k=2;
				for (int i=0;i<lambda_count.size();i++){
					sum_pow+=Math.pow(lambda_count.get(i)[0],k);
				}
			}
		
			

			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			{
				int k=2;
				double sum_local=0;
				for (int i=0;i<lambda_count.size();i++){
					sum_local+=Math.pow(lambda_count.get(i)[0],k);
					if(sum_local/sum_pow>0.99){
						
							break;
					}
					fraction_high+=lambda_count.get(i)[0]/sum_lambda;
					lambda_high.add(lambda_count.get(i)[0]);
				}
			}
		
			double[] lambda_dist=prob_lambda(lambda_count);//0- argument only
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			
			double[] cum_all=new double[lambda_count.size()];
			pp=0;
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]=pp;
				pp+=lambda_count.get(i)[0];
				
			}
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]/=pp;
			}
		
			ArrayList<Double> s1=new ArrayList<Double>();
			ArrayList<Double> s3=new ArrayList<Double>();
			
			
			BinomialDistribution binom=new BinomialDistribution(NN,fraction_high);
			//System.out.println("SART");
			for (int k1=0;k1<10000;k1++){//10000
				int N_high=binom.sample();
				int[] meta_random=random_collision(N_high,cum_high);
				meta_random[0]+=NN-N_high;
				double sum1=0;
				for (int i=1;i<meta_random.length;i++){
					if(meta_random[i]>0){
						sum1+=log_gamma_int(i+1+1)*meta_random[i];
					}
				}
				s1.add(-sum1);
				//for (int i=0;i<meta_random.length;i++){
				//	System.out.print("	"+meta_random[i]);
				//}
				//System.out.println();
			}
			//System.out.println("END");
			
		
			
			
			for (int k2=0;k2<100000;k2++){//100000
				
				s3.add(offset+random_draw(NN,lambda_dist,err,lambda_count));//0 argument only
			//	System.out.println(s3.get(s3.size()-1));
			}
			//System.out.println("END2");
			Collections.sort(s1);
			Collections.sort(s3);
			//System.out.println("END3");
			//System.out.println("Compare");
			/*double sum1=0;
			double sum2=0;
			int nn=0;
			for (int i=0;i<s1.size();i++){
				for (int j=0;j<s3.size();j++){
					if(s1.get(i)+s3.get(j)<=T){
						sum1+=Math.exp(s1.get(i)+s3.get(j));
						nn++;
					}
					sum2+=Math.exp(s1.get(i)+s3.get(j));
				///	System.out.println(i+"	"+j+"	"+Math.exp(s1.get(i)+s3.get(j)));
				}
			}*/
			
			
				int nn=0;
				for (int i=0;i<s3.size();i++){
					if(s3.get(i)+s1.get(0)>T){
						break;
					}
					for (int j=0;j<s1.size();j++){
						if(s3.get(i)+s1.get(j)>T){
							break;
						}
						nn++;
					}
				}
				return (double)(nn)/(double)(s1.size()*s3.size());
			
			
			
			//System.out.println("frac: "+(double)(nn)/(double)(s1.size()*s3.size()));
			
			
			//System.out.println("sum1: "+sum1);
			//System.out.println("sum2: "+sum2);
			//System.out.println("frac_sum: "+(sum1/sum2));
		/*System.out.println("SUM1");
			for (int i=0;i<s1.size();i++){
				System.out.println(Math.exp(s1.get(i)));
			}
			System.out.println("SUM2");
			for (int i=0;i<s3.size();i++){
				System.out.println(Math.exp(s3.get(i)));
			}
			*/
			/*
			System.out.println("Experimental");
			double[] cum=new double[lambda2.size()];
			double p=0;
			for (int i=0;i<cum.length;i++){
				cum[i]=p;
				p+=lambda2.get(i);
			}
			for (int i=0;i<cum.length;i++){
				cum[i]/=p;
			}
			for (int k=0;k<100;k++){
				int[] count=new int[cum.length];
				for (int i=0;i<NN;i++){
					double r=Math.random();
					for (int j=0;j<cum.length;j++){
						if(j<cum.length-1){
							if(cum[j]<=r&&r<cum[j+1]){
								count[j]++;
								break;
							}
						}
						else{
							if(cum[cum.length-1]<=r){
								count[cum.length-1]++;
								break;
							}
						}

					}
				}
				double ss=0;
				for (int i=0;i<count.length;i++){
					if(count[i]>0){
						ss+=Math.log(lambda2.get(i))*count[i];
					}
					if(count[i]>1){
						ss-=Gamma.logGamma(1+count[i]);
					}
					
				}
				System.out.println(Math.exp(ss));
			}
			*/
			
			
			
			
			
			
			
			
			
			
			/*
			double sum1=0;
			double sum2=0;
			double sum3=0;
			for (int i=0;i<s1.size();i++){
				sum1+=(s1.get(i));
			}
			for (int i=0;i<s3.size();i++){
				sum2+=(s3.get(i));
			}
			for (int i=0;i<s3.size();i++){
				if(s3.get(i)+s1.get(0)>T){
					break;
				}
				for (int j=0;j<s1.size();j++){
					if(s3.get(i)+s1.get(j)>T){
						break;
					}
				//	System.out.println(s3.get(i)+s1.get(j));
					sum3+=(s3.get(i)+s1.get(j));
				}
			}*/
			
			/*
			System.out.println("Sum1: "+sum1);
			System.out.println("Sum2: "+sum2);
			System.out.println("Sum3: "+sum3);
			
			return sum3-(s3.size()*sum1+s1.size()*sum2);
			*/
		}
		
		
		public static int[] random_collision (int N, double[] cum){
			int[] count=new int[cum.length];
			for (int i=0;i<N;i++){
				
				double r=random_conc.nextDouble();//Math.rrandom();
				for (int j=0;j<cum.length;j++){
					if(j<cum.length-1){
						if(cum[j]<=r&&r<cum[j+1]){
							count[j]++;
							break;
						}
					}
					else{
						if(cum[cum.length-1]<=r){
							count[cum.length-1]++;
							break;
						}
					}

				}
			}
		
			int[] random_meta=new int[10];
			for (int i=0;i<count.length;i++){
				if(count[i]>0&&count[i]<=random_meta.length){
					random_meta[count[i]-1]++;
				}
			}
			return random_meta;
		}
		
		public static double sign_quick(ArrayList<double[]> lambda_count, int err){
			
			
			int NN=0;
			for (int i=0;i<lambda_count.size();i++){
				NN+=lambda_count.get(i)[1];
			}
			if(lambda_count.size()==0||NN==0){
				return 1;
			}
			
			double T=0;
			for (int i=0;i<lambda_count.size();i++){
				if(lambda_count.get(i)[1]>0){
					T+=lambda_count.get(i)[1]*Math.log(lambda_count.get(i)[0]);
					//System.out.println(Math.log(lambda2.get(i)));
				}
				if(lambda_count.get(i)[1]>0){
					T-=log_gamma_int((int)(lambda_count.get(i)[1])+1);
				}
			}
			double offset=NN*Math.log((double)(1)/(double)(lambda_count.size()));
			T+=offset;
			
			//Comparator<Double> comp=(Double a, Double b)->{
			//	return -a.compareTo(b);
			//};
			//int NN=500;
			
			Collections.sort(lambda_count,comp0);//sort to 0 accordint to count
			double sum_lambda=0;
			for (int i=0;i<lambda_count.size();i++){
				sum_lambda+=lambda_count.get(i)[0];
			}
			
			double sum_pow=0;
			{
				int k=2;
				for (int i=0;i<lambda_count.size();i++){
					sum_pow+=Math.pow(lambda_count.get(i)[0],k);
				}
			}
		
			

			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			{
				int k=2;
				double sum_local=0;
				for (int i=0;i<lambda_count.size();i++){
					sum_local+=Math.pow(lambda_count.get(i)[0],k);
					if(sum_local/sum_pow>0.99){
							break;
					}
					fraction_high+=lambda_count.get(i)[0]/sum_lambda;
					lambda_high.add(lambda_count.get(i)[0]);
				}
			}
		
			double[] lambda_dist=prob_lambda(lambda_count);//according to 0
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			
			double[] cum_all=new double[lambda_count.size()];
			pp=0;
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]=pp;
				pp+=lambda_count.get(i)[0];
				
			}
			for (int i=0;i<cum_all.length;i++){
				cum_all[i]/=pp;
			}
		
			ArrayList<Double> s1=new ArrayList<Double>();
			ArrayList<Double> s3=new ArrayList<Double>();
			
			
			BinomialDistribution binom=new BinomialDistribution(NN,fraction_high);
			//System.out.println("SART");
			for (int k1=0;k1<100;k1++){//
				int N_high=binom.sample();
				int[] meta_random=random_collision(N_high,cum_high);
				meta_random[0]+=NN-N_high;
				double sum1=0;
				for (int i=1;i<meta_random.length;i++){
					if(meta_random[i]>0){
						sum1+=log_gamma_int(i+1+1)*meta_random[i];
					}
				}
				s1.add(-sum1);
				//for (int i=0;i<meta_random.length;i++){
				//	System.out.print("	"+meta_random[i]);
				//}
				//System.out.println();
			}
			//System.out.println("END");
			for (int k2=0;k2<1000;k2++){//
				
				s3.add(offset+random_draw(NN,lambda_dist,err,lambda_count));//0 column
			}
			
			Collections.sort(s1);
			Collections.sort(s3);
			
			int n1=0;
			int n2=0;
			for (int k=0;k<1000;k++){
				if(s3.get((int)(s3.size()*random_conc.nextDouble()))+s1.get((int)(s1.size()*random_conc.nextDouble()))<=T){
					n1++;
				}
				n2++;
			}
			
			if((double)(n1)/(double)(n2)<0.05&&(double)(n1)/(double)(n2)>0.001){
				for (int k=0;k<10000;k++){
					if(s3.get((int)(s3.size()*random_conc.nextDouble()))+s1.get((int)(s1.size()*random_conc.nextDouble()))<=T){
						n1++;
					}
					n2++;
				}
			}
			if((double)(n1)/(double)(n2)<0.005&&(double)(n1)/(double)(n2)>0.001){
				for (int k=0;k<100000;k++){
					if(s3.get((int)(s3.size()*random_conc.nextDouble()))+s1.get((int)(s1.size()*random_conc.nextDouble()))<=T){
						n1++;
					}
					n2++;
				}
			}
					
			
			return (double)(n1)/(double)(n2);
				
			
		}
		
		
		
		static Comparator<double[]> comppp=(double[] a, double[] b)->{
			return -new Double(a[0]).compareTo(new Double(b[0]));
		};
		
		public static ArrayList<double[]> sort (ArrayList<Double> x, ArrayList<Integer> y){
			
			ArrayList<double[]> z=new ArrayList<double[]>();
			for (int i=0;i<x.size();i++){
				z.add(new double[]{x.get(i),y.get(i)});
			}
			Collections.sort(z,comppp);
			return z;
		}
		
		
		public static double hotspot_sign(ArrayList<double[]> pos){
			int max_count=0;
			double lambda_sum=0;
			int count_sum=0;
			
			for (int i=0;i<pos.size();i++){
				if(pos.get(i)[0]>max_count){
					max_count=(int)(pos.get(i)[0]);
				}
				count_sum+=(int)(pos.get(i)[0]);
				lambda_sum+=pos.get(i)[1];
			}
			
			if(max_count>=4&&(double)(max_count)/(double)(count_sum)>=0.05){
				double p=1;
				for (int i=0;i<pos.size();i++){
					p*=new BinomialDistribution(count_sum,pos.get(i)[1]/lambda_sum).cumulativeProbability(max_count-1);
				}
				return 1-p;
			}
			else{
				return 1;
			}
		}
		
		//public static double signREVISED(int count, ArrayList<double[]> lambda2, ArrayList<Integer> count2 ,ArrayList<Double> probability, int err){//int err, 
			//double sign_quick=sign_quickREVISED(count, lambda2,count2, probability,err);//err,
			//if(true){
			//	return sign_quick;
			//}
			//System.out.println(sign_quick);
			/*
			if(sign_quick<0.001){
				double d=sign_longREVISED(count, lambda2,count2,probability,err);
				//if(d<0.00001){
				//	d=sign_longREVISED2(count, lambda2,count2,probability,err);
				//}
				//System.out.println("sign	"+sign_quick+"	"+d);
				return d;//err
			}
			else{
				return sign_quick;
			}*/
			//return sign_quick;
		//}
		
		static int[][] count_quick=new int[100][20000];
		static double[][] s1_quick=new double[100][100]; 
		static double[][] s3_quick=new double[100][1000]; 
		
		static int[][] count_long1=new int[1000][20000];
		static double[][] s1_long1=new double[100][1000]; 
		static double[][] s3_long1=new double[100][10000];
		
		static int[][] count_long2=new int[10000][20000];
		static double[][] s1_long2=new double[100][10000]; 
		static double[][] s3_long2=new double[100][100000];
		
		
		public static double signREVISED(int count_ext, ArrayList<double[]> lambda_count, ArrayList<Double> probability, int err){//int err, 
			//TODO:err
			if(lambda_count.size()==0||count_ext==0||probability.size()==0){
				return 1;
			}
			if(count_ext>=probability.size()||probability.get(count_ext)==0){
				return 0;
			}
			//System.out.println("STARTING Sign Revised");
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<lambda_count.size();i++){
				NN+=lambda_count.get(i)[1];
				sum_lambda+=lambda_count.get(i)[0];
			}
			//System.out.println("NN: "+NN);
			
			if(probability.size()<=NN){
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=log_gamma_int(NN+lambda_count.size()+1);
		
			for (int i=0;i<lambda_count.size();i++){
				if(lambda_count.get(i)[1]>0){
					T+=lambda_count.get(i)[1]*Math.log(lambda_count.get(i)[0]);
					//T1+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(lambda_count.get(i)[1]>1){
					T-=log_gamma_int((int)(lambda_count.get(i)[1])+1);
					//T2+=Gamma.logGamma(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			T+=Math.log(probability.get(NN));
			
			//Comparator<Double> comp=(Double a, Double b)->{
			//	return -a.compareTo(b);
			//};
			
			Collections.sort(lambda_count,comp0);//sort by colum 0 accoridng to comp
			//weitermachen hier
		
			double sum_sq=0;
			for (int i=0;i<lambda_count.size();i++){
				sum_sq+=Math.pow(lambda_count.get(i)[0],2);
			}
		
			double fraction_high=0;
	//		ArrayList<Double> lambda_high=new ArrayList<Double>();
			int lambda_high_n=0;
			double sum_local=0;
			for (int i=0;i<lambda_count.size();i++){
				sum_local+=Math.pow(lambda_count.get(i)[0],2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda_count.get(i)[0]/sum_lambda;
				//lambda_high.add(lambda2.get(i));
				lambda_high_n++;
				if(lambda_high_n>50000){
					return 1;
				}
			}
			
			
			double[] lambda_dist=prob_lambda(lambda_count);//nach 0 colum
			
			
			double[] cum_high=new double[lambda_high_n];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_count.get(i)[0];
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			
			//DDDDDDDDDDDDDD
			
			int[] index_cum_high=create_index(cum_high);
			
			if(probability.size()-1>100){
				//System.out.println("TOO small");
				s1_quick=new double[probability.size()-1][100]; 
				s3_quick=new double[probability.size()-1][1000]; 
				
			}
			
	//		double[][] s1=new double[probability.size()-1][100]; 
	//		double[][] s3=new double[probability.size()-1][1000]; 
			int s1_quick_row=probability.size()-1;
			int s3_quick_row=probability.size()-1;
			int s1_quick_col=100;
			int s3_quick_col=1000;
			
			//System.out.println(lambda_high_n);
			
			//int[][] count_quick=new int[s1_quick_col][lambda_high_n];
			if(lambda_high_n>20000){
				count_quick=new int[s1_quick_col][lambda_high_n];
			}
			for (int i=0;i<s1_quick_col;i++){
				for (int j=0;j<lambda_high_n;j++){
					count_quick[i][j]=0;
				}
			}
			
			//DDDDDDDDDDDDDD
			
			for (int i=0;i<s1_quick_col;i++){
				for (int j=0;j<s1_quick_row;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1_quick[j][i]=0;
						}
						else{
							s1_quick[j][i]=s1_quick[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high,index_cum_high);
						if(count_quick[i][r]==0){
							if(j==0){
								s1_quick[j][i]=0;
							}
							else{
								s1_quick[j][i]=s1_quick[j-1][i];
							}
						}
						else if(count_quick[i][r]==1){
							s1_quick[j][i]=s1_quick[j-1][i]+log_gamma_int(2+1);
						}
						else{
							s1_quick[j][i]=s1_quick[j-1][i]+log_gamma_int(count_quick[i][r]+1+1)-log_gamma_int(count_quick[i][r]+1);
						}
						count_quick[i][r]++;
					}
				}
			}
			
			for (int i=0;i<s3_quick_col;i++){
				for (int j=0;j<s3_quick_row;j++){
					if(j==0){
						s3_quick[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3_quick[j][i]=s3_quick[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
			
			//DDDDDDDDDDDDDDD
			
			/*
			for (int i=0;i<s3_quick_row;i++){
				for (int j=0;j<s3_quick_col;j++){
					for (int k=0;k<err;k++){
						s3_quick[i][j]+=Math.log(lambda_count.get((int)(random_conc.nextDouble()*lambda_count.size()))[0]);
					}
				}
			}*/
			
			//DDDDDDDDDD
			//if(true){
			//	return 1.0;
			//}
			
			/*
			for (int i=0;i<s1_quick_row;i++){
				Arrays.sort(s1_quick[i]);
			}
			for (int i=0;i<s3_quick_row;i++){
				Arrays.sort(s3_quick[i]);
			}*/
			
			////DDDDDDDDDD
			
			double ppp=0;
			double [] cum_probability=new double[probability.size()];
			double sum_prob=0;
			for (int i=0;i<probability.size();i++){
				sum_prob+=probability.get(i);
			}
			for (int i=0;i<probability.size();i++){
				cum_probability[i]=ppp/sum_prob;
				ppp+=probability.get(i);
			}
			int[] index_cum_probability=create_index(cum_probability);
			

			int n1=0;
			int n2=0;
			
			double frac_part=0;
			for (int i=0;i<Math.min(count_ext,probability.size());i++){
				frac_part+=probability.get(i);
			}

			
			for (int k=0;k<1000;k++){
				
				int M=random_draw(cum_probability,count_ext,index_cum_probability);
				
				double exp=log_gamma_int(M+lambda_count.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
				if(M>0){
					exp+=-s1_quick[M-1][(int)(random_conc.nextDouble()*s1_quick_col)]+s3_quick[M-1][(int)(random_conc.nextDouble()*s3_quick_col)]+err(lambda_count,err);
				}
				if(exp<=T){
					n1++;
				}
				n2++;
			}
			
			
			if((double)(n1)/(double)(n2)*(1-frac_part)<0.05&&(double)(n1)/(double)(n2)*(1-frac_part)>0.001){
				for (int k=0;k<10000;k++){
					
					int M=random_draw(cum_probability,count_ext,index_cum_probability);
					
					double exp=log_gamma_int(M+lambda_count.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
					if(M>0){
						exp+=-s1_quick[M-1][(int)(random_conc.nextDouble()*s1_quick_col)]+s3_quick[M-1][(int)(random_conc.nextDouble()*s3_quick_col)]+err(lambda_count,err);
					}
					if(exp<=T){
						n1++;
					}
					n2++;
				}
			}
			
			if((double)(n1)/(double)(n2)*(1-frac_part)<0.005&&(double)(n1)/(double)(n2)*(1-frac_part)>0.001){
				for (int k=0;k<100000;k++){
					
					int M=random_draw(cum_probability,count_ext,index_cum_probability);
					
					double exp=log_gamma_int(M+lambda_count.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
					if(M>0){
						exp+=-s1_quick[M-1][(int)(random_conc.nextDouble()*s1_quick_col)]+s3_quick[M-1][(int)(random_conc.nextDouble()*s3_quick_col)]+err(lambda_count,err);
					}
					if(exp<=T){
						n1++;
					}
					n2++;
				}
			}
			
			if(s1_quick.length>100){
				s1_quick=new double[100][100]; 
				s3_quick=new double[100][1000];
			}
			
			if(count_quick[0].length>20000){
				count_quick=new int[s1_quick_col][20000];
			}
			
			double sign= (double)(n1)/(double)(n2)*(1-frac_part);
			
			if(sign>=0.001){
				return sign;
			}
			
			
			sign=0;
			{
				if(probability.size()-1>s1_long1.length){
					s1_long1=new double[probability.size()-1][1000]; //10000
					s3_long1=new double[probability.size()-1][10000];//100000
				}
				
				//double[][] s1=new double[probability.size()-1][10000]; 
				//double[][] s3=new double[probability.size()-1][100000];
				int s1_col=1000;//10000
				int s3_col=10000;//100000
				int s1_row=probability.size()-1;
				int s3_row=probability.size()-1;
				
				
				
//				System.out.println("generating arrays done");
				
				//System.out.println(lambda_high.size());//TODO: 
				if(lambda_high_n>20000){
					count_long1=new int[s1_col][lambda_high_n];
				}
				
				for (int i=0;i<s1_col;i++){
					for (int j=0;j<lambda_high_n;j++){
						count_long1[i][j]=0;
					}
				}
				
				//count_long1=new int[s1_col][lambda_high_n];
				for (int i=0;i<s1_col;i++){
					for (int j=0;j<s1_row;j++){
						if(random_conc.nextDouble()>fraction_high){
							if(j==0){
								s1_long1[j][i]=0;
							}
							else{
								s1_long1[j][i]=s1_long1[j-1][i];
							}
						}
						else{
							int r=random_draw(cum_high,index_cum_high);
							if(count_long1[i][r]==0){
								if(j==0){
									s1_long1[j][i]=0;
								}
								else{
									s1_long1[j][i]=s1_long1[j-1][i];
								}
							}
							else if(count_long1[i][r]==1){
								s1_long1[j][i]=s1_long1[j-1][i]-log_gamma_int(2+1);
							}
							else{
								s1_long1[j][i]=s1_long1[j-1][i]-log_gamma_int(count_long1[i][r]+1+1)+log_gamma_int(count_long1[i][r]+1);
							}
							count_long1[i][r]++;
						}
					}
				}
//				System.out.println("s1 done");
				
				for (int i=0;i<s3_col;i++){
					for (int j=0;j<s3_row;j++){
						if(j==0){
							s3_long1[j][i]=random_draw_log(lambda_dist);
						}
						else{
							s3_long1[j][i]=s3_long1[j-1][i]+random_draw_log(lambda_dist);
						}
					}
				}
//				System.out.println("s3 done");
				for (int i=0;i<s3_row;i++){
					for (int j=0;j<s3_col;j++){
						for (int k=0;k<err;k++){
							s3_long1[i][j]+=Math.log(lambda_count.get((int)(random_conc.nextDouble()*lambda_count.size()))[0]);
						}
					}
				}
				
				
				for (int i=0;i<s1_row;i++){
					Arrays.sort(s1_long1[i]);
				}
				for (int i=0;i<s3_row;i++){
					Arrays.sort(s3_long1[i]);
				}
				
				
				for (int M=count_ext;M<probability.size();M++){
					if(probability.get(M)==0){
						continue;
					}
					double start=+log_gamma_int(M+lambda_count.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
					double frac=1;
					if(M>0){
						int nn=0;
						for (int i=0;i<1000;i++){
							if(s3_long1[M-1][(int)(random_conc.nextDouble()*s3_col)]+s1_long1[M-1][(int)(s1_col*random_conc.nextDouble())]+start<=T){
								nn++;
							}
						}
						double f=(double)(nn)/1000.0;
						if(f<=0.05){
							
							//System.out.println("Intense "+f);
							nn=0;
							for (int i=0;i<s3_col;i++){
								if(s3_long1[M-1][i]+s1_long1[M-1][0]+start>T){
									break;
								}
								for (int j=0;j<s1_col;j++){
									if(s3_long1[M-1][i]+s1_long1[M-1][j]+start>T){
										break;
									}
									nn++;
								}
							}
							frac=(double)(nn)/(double)(s1_col*s3_col);
							
						
						}
						else{
							frac=f;
						}
						
					}
					else{
						if(Math.log(probability.get(M))<=T){
							frac=1;
						}
						else{
							frac=0;
						}
					}
					sign+=probability.get(M)*frac;
					
				}
				
				if(s1_long1.length>100){
					s1_long1=new double[100][1000]; //10000
					s3_long1=new double[100][10000];//100000
					
				}
				if(count_long1[0].length>20000){
					count_long1=new int[s1_col][20000];
				}
				
			}
			
			if(sign<0.00001){
				sign=0;
				
				if(probability.size()-1>s1_long2.length){
					s1_long2=new double[probability.size()-1][10000]; 
					s3_long2=new double[probability.size()-1][100000];
				}
				
				int s1_col=10000;
				int s3_col=100000;
				int s1_row=probability.size()-1;
				int s3_row=probability.size()-1;
				
				
				
//				System.out.println("generating arrays done");
				
				//System.out.println(lambda_high.size());//TODO: 
				
				if(lambda_high_n>20000){
					count_long2=new int[s1_col][lambda_high_n];
				}
				
				for (int i=0;i<s1_col;i++){
					for (int j=0;j<lambda_high_n;j++){
						count_long2[i][j]=0;
					}
				}
				
				//count=new int[s1_col][lambda_high_n];
				for (int i=0;i<s1_col;i++){
					for (int j=0;j<s1_row;j++){
						if(random_conc.nextDouble()>fraction_high){
							if(j==0){
								s1_long2[j][i]=0;
							}
							else{
								s1_long2[j][i]=s1_long2[j-1][i];
							}
						}
						else{
							int r=random_draw(cum_high,index_cum_high);
							if(count_long2[i][r]==0){
								if(j==0){
									s1_long2[j][i]=0;
								}
								else{
									s1_long2[j][i]=s1_long2[j-1][i];
								}
							}
							else if(count_long2[i][r]==1){
								s1_long2[j][i]=s1_long2[j-1][i]-log_gamma_int(2+1);
							}
							else{
								s1_long2[j][i]=s1_long2[j-1][i]-log_gamma_int(count_long2[i][r]+1+1)+log_gamma_int(count_long2[i][r]+1);
							}
							count_long2[i][r]++;
						}
					}
				}
//				System.out.println("s1 done");
				
				for (int i=0;i<s3_col;i++){
					for (int j=0;j<s3_row;j++){
						if(j==0){
							s3_long2[j][i]=random_draw_log(lambda_dist);
						}
						else{
							s3_long2[j][i]=s3_long2[j-1][i]+random_draw_log(lambda_dist);
						}
					}
				}
//				System.out.println("s3 done");
				for (int i=0;i<s3_row;i++){
					for (int j=0;j<s3_col;j++){
						for (int k=0;k<err;k++){
							s3_long2[i][j]+=Math.log(lambda_count.get((int)(random_conc.nextDouble()*lambda_count.size()))[0]);
						}
					}
				}
				
				
				for (int i=0;i<s1_row;i++){
					Arrays.sort(s1_long2[i]);
				}
				for (int i=0;i<s3_row;i++){
					Arrays.sort(s3_long2[i]);
				}
				
				
				for (int M=count_ext;M<probability.size();M++){
					if(probability.get(M)==0){
						continue;
					}
					double start=+log_gamma_int(M+lambda_count.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
					double frac=1;
					if(M>0){
						int nn=0;
						for (int i=0;i<1000;i++){
							if(s3_long2[M-1][(int)(random_conc.nextDouble()*s3_col)]+s1_long2[M-1][(int)(s1_col*random_conc.nextDouble())]+start<=T){
								nn++;
							}
						}
						double f=(double)(nn)/1000.0;
						if(f<=0.05){
							
							//System.out.println("Intense "+f);
							nn=0;
							for (int i=0;i<s3_col;i++){
								if(s3_long2[M-1][i]+s1_long2[M-1][0]+start>T){
									break;
								}
								for (int j=0;j<s1_col;j++){
									if(s3_long2[M-1][i]+s1_long2[M-1][j]+start>T){
										break;
									}
									nn++;
								}
							}
							frac=(double)(nn)/(double)(s1_col*s3_col);
							
						
						}
						else{
							frac=f;
						}
						
					}
					else{
						if(Math.log(probability.get(M))<=T){
							frac=1;
						}
						else{
							frac=0;
						}
					}
					sign+=probability.get(M)*frac;
					
				}
				
				if(s1_long2.length>100){
					s1_long2=new double[100][10000]; 
					s3_long2=new double[100][100000];
				}
				if(count_long2[0].length>20000){
					count_long2=new int[s1_col][20000];
				}
				
			}
			

			return sign;
			
		}
		
		public static double err(ArrayList<double[]> lambda_count, int err){
			double sum=0;
			for (int k=0;k<err;k++){
				sum+=Math.log(lambda_count.get((int)(random_conc.nextDouble()*lambda_count.size()))[0]);
			}
			return sum;
		}
		public static double[] prob_lambda(ArrayList<double[]> lambda_high){
			if(lambda_high.size()==0){
				return new double[0];
			}
			double sum=0;
			for (int i=0;i<lambda_high.size();i++){
				sum+=lambda_high.get(i)[0];
			}
			
			double[] array=new double[100000];
			for (int i=0;i<array.length;i++){
				array[i]=Double.NaN;
			}
			double cum=0;
			for (int i=0;i<lambda_high.size();i++){
				array[(int)(cum*(array.length-1))]=Math.log(lambda_high.get(i)[0]);
				cum+=lambda_high.get(i)[0]/sum;
			}
			
			double prev=Double.NaN;
			for (int i=array.length-1;i>=0;i--){
				if(Double.isNaN(array[i])){
					array[i]=prev;
				}
				else{
					prev=array[i];
				}
			}
			int i=array.length-1;
			while(i>=0&&Double.isNaN(array[i])){
				i--;
			}
			for (int ii=i+1;ii<array.length;ii++){
				array[ii]=array[i];
			}
			
			return array;
			
			
		}
		
		public static int[] create_index(double[] cum){
			int[] index=new int[101];
			for (int k=0;k<=100;k++){
				double r=(double)(k)/100.0;
				for (int i=0;i<cum.length;i++){
					if(i<cum.length-1){
						if(cum[i]<=r&&r<cum[i+1]){
							index[k]=i;
							break;
						}
					}
					else{
						if(cum[i]<=r){
							index[k]=i;
							break;
						}
					}
				}
			}
			//for (int k=0;k<index.length;k++){
			//	System.out.println(k+"	"+index[k]+"	"+cum[index[k]]);
			//}
			
			
			return index;
		}
		
		/*
		public static double sign_longREVISED(int count_ext, ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability, int err){//int err, 
			//TODO err
			if(lambda2.size()==0||count_ext==0||probability.size()==0){
				return 1;
			}
			if(count_ext>=probability.size()||probability.get(count_ext)==0){
				return 0;
			}
			
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				System.out.println("Probability too short");
				System.exit(0);
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=log_gamma_int(NN+count2.size()+1);
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=log_gamma_int(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			T+=Math.log(probability.get(NN));
			
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
				if(lambda_high.size()>50000){
					return 1;
				}
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			int[] index_cum_high=create_index(cum_high);
			
//			System.out.println("generating arrays");
			
			double[][] s1=new double[probability.size()-1][10000]; 
			double[][] s3=new double[probability.size()-1][100000]; 
			
//			System.out.println("generating arrays done");
			
			//System.out.println(lambda_high.size());//TODO: 
			int[][] count=new int[s1[0].length][lambda_high.size()];
			for (int i=0;i<s1[0].length;i++){
				for (int j=0;j<s1.length;j++){
					if(random_conc.nextDouble()>fraction_high){
						if(j==0){
							s1[j][i]=0;
						}
						else{
							s1[j][i]=s1[j-1][i];
						}
					}
					else{
						int r=random_draw(cum_high,index_cum_high);
						if(count[i][r]==0){
							if(j==0){
								s1[j][i]=0;
							}
							else{
								s1[j][i]=s1[j-1][i];
							}
						}
						else if(count[i][r]==1){
							s1[j][i]=s1[j-1][i]-log_gamma_int(2+1);
						}
						else{
							s1[j][i]=s1[j-1][i]-log_gamma_int(count[i][r]+1+1)+log_gamma_int(count[i][r]+1);
						}
						count[i][r]++;
					}
				}
			}
//			System.out.println("s1 done");
			
			for (int i=0;i<s3[0].length;i++){
				for (int j=0;j<s3.length;j++){
					if(j==0){
						s3[j][i]=random_draw_log(lambda_dist);
					}
					else{
						s3[j][i]=s3[j-1][i]+random_draw_log(lambda_dist);
					}
				}
			}
//			System.out.println("s3 done");
			for (int i=0;i<s3.length;i++){
				for (int j=0;j<s3[i].length;j++){
					for (int k=0;k<err;k++){
						s3[i][j]+=Math.log(lambda2.get((int)(random_conc.nextDouble()*lambda2.size())));
					}
				}
			}
			
			
			for (int i=0;i<s1.length;i++){
				Arrays.sort(s1[i]);
			}
			for (int i=0;i<s3.length;i++){
				Arrays.sort(s3[i]);
			}
			
			double sign=0;
			
			for (int M=count_ext;M<probability.size();M++){
				if(probability.get(M)==0){
					continue;
				}
				double start=+log_gamma_int(M+count2.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
				double frac=1;
				if(M>0){
					int nn=0;
					for (int i=0;i<1000;i++){
						if(s3[M-1][(int)(random_conc.nextDouble()*s3[M-1].length)]+s1[M-1][(int)(s1[M-1].length*random_conc.nextDouble())]+start<=T){
							nn++;
						}
					}
					double f=(double)(nn)/1000.0;
					if(f<=0.05){
						
						//System.out.println("Intense "+f);
						nn=0;
						for (int i=0;i<s3[M-1].length;i++){
							if(s3[M-1][i]+s1[M-1][0]+start>T){
								break;
							}
							for (int j=0;j<s1[M-1].length;j++){
								if(s3[M-1][i]+s1[M-1][j]+start>T){
									break;
								}
								nn++;
							}
						}
						frac=(double)(nn)/(double)(s1[M-1].length*s3[M-1].length);
						
					
					}
					else{
						frac=f;
					}
					
				}
				else{
					if(Math.log(probability.get(M))<=T){
						frac=1;
					}
					else{
						frac=0;
					}
				}
				sign+=probability.get(M)*frac;
				
			}
			
			
			return sign;

			
		}*/
		
		
		/*
		public static double sign_longREVISED(int count_ext, ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability, int err){//int err, 
			//TODO err
			if(lambda2.size()==0||count_ext==0||probability.size()==0){
				return 1;
			}
			if(count_ext>=probability.size()||probability.get(count_ext)==0){
				return 0;
			}
			
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				System.out.println("Probability too short");
				System.exit(0);
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=log_gamma_int(NN+count2.size()+1);
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=log_gamma_int(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			T+=Math.log(probability.get(NN));
			
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
				if(lambda_high.size()>50000){
					return 1;
				}
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			int[] index_cum_high=create_index(cum_high);
			
//			System.out.println("generating arrays");
			
			double sign=0;
			{
				if(probability.size()-1>s1_long1.length){
					s1_long1=new double[probability.size()-1][1000]; //10000
					s3_long1=new double[probability.size()-1][10000];//100000
				}
				
				//double[][] s1=new double[probability.size()-1][10000]; 
				//double[][] s3=new double[probability.size()-1][100000];
				int s1_col=1000;//10000
				int s3_col=10000;//100000
				int s1_row=probability.size()-1;
				int s3_row=probability.size()-1;
				
				
				
//				System.out.println("generating arrays done");
				
				//System.out.println(lambda_high.size());//TODO: 
				int[][] count=new int[s1_col][lambda_high.size()];
				for (int i=0;i<s1_col;i++){
					for (int j=0;j<s1_row;j++){
						if(random_conc.nextDouble()>fraction_high){
							if(j==0){
								s1_long1[j][i]=0;
							}
							else{
								s1_long1[j][i]=s1_long1[j-1][i];
							}
						}
						else{
							int r=random_draw(cum_high,index_cum_high);
							if(count[i][r]==0){
								if(j==0){
									s1_long1[j][i]=0;
								}
								else{
									s1_long1[j][i]=s1_long1[j-1][i];
								}
							}
							else if(count[i][r]==1){
								s1_long1[j][i]=s1_long1[j-1][i]-log_gamma_int(2+1);
							}
							else{
								s1_long1[j][i]=s1_long1[j-1][i]-log_gamma_int(count[i][r]+1+1)+log_gamma_int(count[i][r]+1);
							}
							count[i][r]++;
						}
					}
				}
//				System.out.println("s1 done");
				
				for (int i=0;i<s3_col;i++){
					for (int j=0;j<s3_row;j++){
						if(j==0){
							s3_long1[j][i]=random_draw_log(lambda_dist);
						}
						else{
							s3_long1[j][i]=s3_long1[j-1][i]+random_draw_log(lambda_dist);
						}
					}
				}
//				System.out.println("s3 done");
				for (int i=0;i<s3_row;i++){
					for (int j=0;j<s3_col;j++){
						for (int k=0;k<err;k++){
							s3_long1[i][j]+=Math.log(lambda2.get((int)(random_conc.nextDouble()*lambda2.size())));
						}
					}
				}
				
				
				for (int i=0;i<s1_row;i++){
					Arrays.sort(s1_long1[i]);
				}
				for (int i=0;i<s3_row;i++){
					Arrays.sort(s3_long1[i]);
				}
				
				
				for (int M=count_ext;M<probability.size();M++){
					if(probability.get(M)==0){
						continue;
					}
					double start=+log_gamma_int(M+count2.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
					double frac=1;
					if(M>0){
						int nn=0;
						for (int i=0;i<1000;i++){
							if(s3_long1[M-1][(int)(random_conc.nextDouble()*s3_long1[M-1].length)]+s1_long1[M-1][(int)(s1_long1[M-1].length*random_conc.nextDouble())]+start<=T){
								nn++;
							}
						}
						double f=(double)(nn)/1000.0;
						if(f<=0.05){
							
							//System.out.println("Intense "+f);
							nn=0;
							for (int i=0;i<s3_long1[M-1].length;i++){
								if(s3_long1[M-1][i]+s1_long1[M-1][0]+start>T){
									break;
								}
								for (int j=0;j<s1_long1[M-1].length;j++){
									if(s3_long1[M-1][i]+s1_long1[M-1][j]+start>T){
										break;
									}
									nn++;
								}
							}
							frac=(double)(nn)/(double)(s1_long1[M-1].length*s3_long1[M-1].length);
							
						
						}
						else{
							frac=f;
						}
						
					}
					else{
						if(Math.log(probability.get(M))<=T){
							frac=1;
						}
						else{
							frac=0;
						}
					}
					sign+=probability.get(M)*frac;
					
				}
				
				if(s1_long1.length>100){
					s1_long1=new double[100][1000]; //10000
					s3_long1=new double[100][10000];//100000
					
				}
				
				
			}
			
			if(sign<0.00001){
				sign=0;
				
				if(probability.size()-1>s1_long2.length){
					s1_long2=new double[probability.size()-1][10000]; 
					s3_long2=new double[probability.size()-1][100000];
				}
				
				int s1_col=10000;
				int s3_col=100000;
				int s1_row=probability.size()-1;
				int s3_row=probability.size()-1;
				
				
				
//				System.out.println("generating arrays done");
				
				//System.out.println(lambda_high.size());//TODO: 
				int[][] count=new int[s1_col][lambda_high.size()];
				for (int i=0;i<s1_col;i++){
					for (int j=0;j<s1_row;j++){
						if(random_conc.nextDouble()>fraction_high){
							if(j==0){
								s1_long2[j][i]=0;
							}
							else{
								s1_long2[j][i]=s1_long2[j-1][i];
							}
						}
						else{
							int r=random_draw(cum_high,index_cum_high);
							if(count[i][r]==0){
								if(j==0){
									s1_long2[j][i]=0;
								}
								else{
									s1_long2[j][i]=s1_long2[j-1][i];
								}
							}
							else if(count[i][r]==1){
								s1_long2[j][i]=s1_long2[j-1][i]-log_gamma_int(2+1);
							}
							else{
								s1_long2[j][i]=s1_long2[j-1][i]-log_gamma_int(count[i][r]+1+1)+log_gamma_int(count[i][r]+1);
							}
							count[i][r]++;
						}
					}
				}
//				System.out.println("s1 done");
				
				for (int i=0;i<s3_col;i++){
					for (int j=0;j<s3_row;j++){
						if(j==0){
							s3_long2[j][i]=random_draw_log(lambda_dist);
						}
						else{
							s3_long2[j][i]=s3_long2[j-1][i]+random_draw_log(lambda_dist);
						}
					}
				}
//				System.out.println("s3 done");
				for (int i=0;i<s3_row;i++){
					for (int j=0;j<s3_col;j++){
						for (int k=0;k<err;k++){
							s3_long2[i][j]+=Math.log(lambda2.get((int)(random_conc.nextDouble()*lambda2.size())));
						}
					}
				}
				
				
				for (int i=0;i<s1_row;i++){
					Arrays.sort(s1_long2[i]);
				}
				for (int i=0;i<s3_row;i++){
					Arrays.sort(s3_long2[i]);
				}
				
				
				for (int M=count_ext;M<probability.size();M++){
					if(probability.get(M)==0){
						continue;
					}
					double start=+log_gamma_int(M+count2.size()+1)-M*Math.log(sum_lambda)+Math.log(probability.get(M));
					double frac=1;
					if(M>0){
						int nn=0;
						for (int i=0;i<1000;i++){
							if(s3_long2[M-1][(int)(random_conc.nextDouble()*s3_long2[M-1].length)]+s1_long2[M-1][(int)(s1_long2[M-1].length*random_conc.nextDouble())]+start<=T){
								nn++;
							}
						}
						double f=(double)(nn)/1000.0;
						if(f<=0.05){
							
							//System.out.println("Intense "+f);
							nn=0;
							for (int i=0;i<s3_long2[M-1].length;i++){
								if(s3_long2[M-1][i]+s1_long2[M-1][0]+start>T){
									break;
								}
								for (int j=0;j<s1_long2[M-1].length;j++){
									if(s3_long2[M-1][i]+s1_long2[M-1][j]+start>T){
										break;
									}
									nn++;
								}
							}
							frac=(double)(nn)/(double)(s1_long2[M-1].length*s3_long2[M-1].length);
							
						
						}
						else{
							frac=f;
						}
						
					}
					else{
						if(Math.log(probability.get(M))<=T){
							frac=1;
						}
						else{
							frac=0;
						}
					}
					sign+=probability.get(M)*frac;
					
				}
				
				if(s1_long2.length>100){
					s1_long2=new double[100][10000]; 
					s3_long2=new double[100][100000];
				}
				
				
			}
			

			return sign;
		}
	*/	
		
	/*	
		public static double sign_longREVISED2(int count_ext, ArrayList<Double> lambda2, ArrayList<Integer> count2, ArrayList<Double> probability, int err){//int err, 
			//TODO err
			if(lambda2.size()==0||count_ext==0||probability.size()==0){
				return 1;
			}
			if(count_ext>=probability.size()||probability.get(count_ext)==0){
				return 0;
			}
			
			if(probability.size()==0){
				return 1;
			}
			int NN=0;
			double sum_lambda=0;
			for (int i=0;i<count2.size();i++){
				NN+=count2.get(i);
				sum_lambda+=lambda2.get(i);
			}
			if(lambda2.size()==0||NN==0){
				return 1;
			}
			if(probability.size()<=NN){
				System.out.println("Probability too short");
				System.exit(0);
				return 0;
			}
			if(probability.get(NN)==0){
				return 0;
			}
			
			double T=log_gamma_int(NN+count2.size()+1);
			for (int i=0;i<count2.size();i++){
				if(count2.get(i)>0){
					T+=count2.get(i)*Math.log(lambda2.get(i));
				}
				if(count2.get(i)>1){
					T-=log_gamma_int(count2.get(i)+1);
				}
			}
			T-=Math.log(sum_lambda)*NN;
			T+=Math.log(probability.get(NN));
			
			
			Comparator<Double> comp=(Double a, Double b)->{
				return -a.compareTo(b);
			};
			
			Collections.sort(lambda2,comp);
			
		
			double sum_sq=0;
			for (int i=0;i<lambda2.size();i++){
				sum_sq+=Math.pow(lambda2.get(i),2);
			}
		
			double fraction_high=0;
			ArrayList<Double> lambda_high=new ArrayList<Double>();
			double sum_local=0;
			for (int i=0;i<lambda2.size();i++){
				sum_local+=Math.pow(lambda2.get(i),2);
				if(sum_local/sum_sq>0.99){
						break;
				}
				fraction_high+=lambda2.get(i)/sum_lambda;
				lambda_high.add(lambda2.get(i));
				if(lambda_high.size()>50000){
					return 1;
				}
			}
			
			double[] lambda_dist=prob_lambda(lambda2);
			
			
			double[] cum_high=new double[lambda_high.size()];
			double pp=0;
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]=pp;
				pp+=lambda_high.get(i);
				
			}
			for (int i=0;i<cum_high.length;i++){
				cum_high[i]/=pp;
			}
			int[] index_cum_high=create_index(cum_high);
			
//			System.out.println("generating arrays");
			
			

			
		}*/
		
		public static double random_draw_log(double[] lambda){
			
			return lambda[(int)(lambda.length*random_conc.nextDouble())];
			
			
		}
		
		public static double random_draw(int N, double[] lambda, int N_err, ArrayList<double[]> lambda2){
			if(N==0){
				return 0;
			}
			double sum=0;
			for (int i=0;i<N;i++){
				sum+=lambda[(int)(lambda.length*random_conc.nextDouble())];
			}
			for (int i=0;i<N_err;i++){
				sum+=Math.log(lambda2.get((int)(random_conc.nextDouble()*lambda2.size()))[0]);
			}
			return sum;
		}
		
		
		public static int random_draw(double[] cum, int min, int[] index){
			
			double r=(1-cum[min])*random_conc.nextDouble()+cum[min];
			for (int i=Math.max(0, index[(int)(r*100)]-1);i<Math.min(cum.length,index[Math.min(index.length-1, (int)(r*100)+1)]+1);i++){
			//for (int i=0;i<cum.length;i++){
				if(i<cum.length-1){
					if(cum[i]<=r&&r<cum[i+1]){
						return i;
					}
				}
				else{
					if(cum[i]<=r){
						return i;
					}
				}
			}
			return -1;
		}
		
		public static int random_draw(double[] cum, int[] index){
			double r=random_conc.nextDouble();
			//for (int i=0;i<cum.length;i++){
			for (int i=Math.max(0,index[(int)(r*100)]-1);i<Math.min(index[(int)(r*100)+1]+1,cum.length);i++){
				if(i<cum.length-1){
					if(cum[i]<=r&&r<cum[i+1]){
						//System.out.println(i+"	vs	"+random_draw(cum,r));
						return i;
					}
				}
				else{
					if(cum[i]<=r){
						//System.out.println(i+"	vs	"+random_draw(cum,r));
						return i;
					}
				}
			}
			return -1;
		}
		
		
		public static int random_draw(double[] cum, int min, double r){
			
			for (int i=0;i<cum.length;i++){
				if(i<cum.length-1){
					if(cum[i]<=r&&r<cum[i+1]){
						//System.out.println(i+"	vs	"+random_draw(cum,min,r));
						return i;
					}
				}
				else{
					if(cum[i]<=r){
						//System.out.println(i+"	vs	"+random_draw(cum,min,r));
						return i;
					}
				}
			}
			return -1;
		}
		
		public static int random_draw(double[] cum, double r){
			for (int i=0;i<cum.length;i++){
				if(i<cum.length-1){
					if(cum[i]<=r&&r<cum[i+1]){
						return i;
					}
				}
				else{
					if(cum[i]<=r){
						return i;
					}
				}
			}
			return -1;
		}
		
		
		public static ArrayList<double[]> summarize(ArrayList<int[][]> count,ArrayList<double[][]> lambda2,ArrayList<String> nucl2, ArrayList<String> amino_acid, ArrayList<int[]> index_gene, int k){
			boolean[] taken=new boolean[index_gene.size()];
			ArrayList<double[]> pos2=new ArrayList<double[]>();
			for (int ii=0;ii<index_gene.size();ii++){
				if(taken[ii]){
					continue;
				}
				if(amino_acid.get(index_gene.get(ii)[0]).equals("-")){
					
					double x=0;
					if(lambda2.get(index_gene.get(ii)[0]).length==1){
						int iii=(int)(lambda2.get(index_gene.get(ii)[0])[0][0]);
						if(nucl2.get(index_gene.get(ii)[0]).equals("C")||nucl2.get(index_gene.get(ii)[0]).equals("G")){
							x=lambda_context_product6_weight[k][tt2[index_gene.get(ii)[1]]][iii];
						}
						else {
							x=lambda_context_product6_weight[k][tt1[index_gene.get(ii)[1]]][iii];
						}
					}
					else{
						x=lambda2.get(index_gene.get(ii)[0])[k][index_gene.get(ii)[1]];
					}
					
					
					pos2.add(new double[]{count.get(index_gene.get(ii)[0])[k][index_gene.get(ii)[1]],x});//index_gene.get(ii)[1]//pos.get(i)
					taken[ii]=true;
					continue;
				}
				
				double[] a=new double[2];
				for (int j=0;j<9;j++){
					if(ii+j>=index_gene.size()){
						break;
					}
					if(!amino_acid.get(index_gene.get(ii+j)[0]).equals("-")&&amino_acid.get(index_gene.get(ii+j)[0]).equals(amino_acid.get(index_gene.get(ii)[0]))){
						a[0]+=count.get(index_gene.get(ii+j)[0])[k][index_gene.get(ii+j)[1]];//pos.get(i+j)[0];
						//a[1]+=lambda.get(index_gene.get(ii+j)[0])[index_gene.get(ii+j)[1]];//pos.get(i+j)[1];
						
						//double x=0;
						if(lambda2.get(index_gene.get(ii+j)[0]).length==1){
							int iii=(int)(lambda2.get(index_gene.get(ii+j)[0])[0][0]);
							if(nucl2.get(index_gene.get(ii+j)[0]).equals("C")||nucl2.get(index_gene.get(ii+j)[0]).equals("G")){
								a[1]+=lambda_context_product6_weight[k][tt2[index_gene.get(ii+j)[1]]][iii];
							}
							else {
								a[1]+=lambda_context_product6_weight[k][tt1[index_gene.get(ii+j)[1]]][iii];
							}
						}
						else{
							a[1]+=lambda2.get(index_gene.get(ii+j)[0])[k][index_gene.get(ii+j)[1]];
						}
						
						
						
						
						taken[ii+j]=true;
					}
				}
				pos2.add(a);
			}
			return pos2;
		
		}
		
		public static ArrayList<double[]> summarize(ArrayList<double[]> pos, ArrayList<String> aa){
			boolean[] taken=new boolean[pos.size()];
			ArrayList<double[]> pos2=new ArrayList<double[]>();
			for (int i=0;i<pos.size();i++){
				if(taken[i]){
					continue;
				}
				if(aa.get(i).equals("-")){
					pos2.add(pos.get(i));
					taken[i]=true;
					continue;
				}
				
				double[] a=new double[2];
				for (int j=0;j<9;j++){
					if(i+j>=pos.size()){
						break;
					}
					if(!aa.get(i+j).equals("-")&&aa.get(i+j).equals(aa.get(i))){
						a[0]+=pos.get(i+j)[0];
						a[1]+=pos.get(i+j)[1];
						taken[i+j]=true;
					}
				}
				pos2.add(a);
			}
			return pos2;
		}
		
		
		
		
	
		
		//the "heart" of this thread which computes the local mutation rate in the center of a queu and appends it to output
		//in brief walks from -10 to +10 around the central position (skipping 0) and forming the product over likelihoods
		
		
		
		
		static int[] tt1={3,5,4};
		static int[] tt2={0,1,2};
		static int offset_left=10;
		static int offset_right=10;
		public static void update_10_10(){
			try{
				//int x=0;
				//int[] xx=new int[20];
				
				
				boolean valid=true;
				int index=0;
				if(nucl.get(50).equals("C")||nucl.get(50).equals("T")){
					int n=1;
					for (int j=-5;j<=-1;j++){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){
							index+=n*nucl_index.get(50+j);
							n*=4;
						}
						else{
							valid=false;
						}
						
					}
					for (int j=1;j<=5;j++){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){
							index+=n*nucl_index.get(50+j);
							n*=4;
						}
						else{
							valid=false;
						}
					}
				}
				else {
					int n=1;
					for (int j=5;j>=1;j--){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){
							index+=n*(3-nucl_index.get(50+j));
							n*=4;
						}
						else{
							valid=false;
						}
						
					}
					for (int j=-1;j>=-5;j--){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){
							index+=n*(3-nucl_index.get(50+j));
							n*=4;
						}
						else{
							valid=false;
						}
					}
				}
				
				pos2.add(pos.get(50));
				nucl2.add(nucl.get(50));
				coverage2.add(coverage.get(50));
				label2.add(label.get(50));
				count2.add(count.get(50));
				amino_acid2.add(amino_acid.get(50));
				
				if (valid){
					lambda2.add(new double[][]{{index}});
				
				}
				else{
					double[][] lambda=new double[lambda_context.size()][3];
					
					if(nucl.get(50).equals("A")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][1];
							lambda[i][1]*=lambda_type.get(i)[0][1];
							lambda[i][2]*=lambda_type.get(i)[0][1];
							for (int j=-1;j>=-5;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){
											lambda[i][k]*=lambda_context.get(i)[9-j][tt1[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							for (int j=1;j<=5;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10-j][tt1[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}		
							}
						}
					}
					else if(nucl.get(50).equals("C")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][0];
							lambda[i][1]*=lambda_type.get(i)[0][0];
							lambda[i][2]*=lambda_type.get(i)[0][0];
							
							for (int j=-1;j>=-5;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10+j][tt2[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							for (int j=1;j<=5;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){
											lambda[i][k]*=lambda_context.get(i)[9+j][tt2[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
						}
					}
					else if(nucl.get(50).equals("G")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][0];
							lambda[i][1]*=lambda_type.get(i)[0][0];
							lambda[i][2]*=lambda_type.get(i)[0][0];
							
							for (int j=-1;j>=-5;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){		
											lambda[i][k]*=lambda_context.get(i)[9-j][tt2[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;	
								}
							}
							for (int j=1;j<=5;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10-j][tt2[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
								
							}
						}
					}
					else if(nucl.get(50).equals("T")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][1];
							lambda[i][1]*=lambda_type.get(i)[0][1];
							lambda[i][2]*=lambda_type.get(i)[0][1];
							
							
							for (int j=-1;j>=-5;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10+j][tt1[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							for (int j=1;j<=5;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){	
											lambda[i][k]*=lambda_context.get(i)[9+j][tt1[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							
						}
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					
					
					//if(!valid){
						lambda2.add(product(weights,weights_index,lambda));
						
					//}
					//lambda22.add(product(weights,weights_index,lambda));
					
					
				}
				
				/*
				for (int j=-1;j>=-offset_left;j--){
					if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){//nucl_index(nucl.get(50+j))!=-1){
						//xx[10+j]=nucl_index(nucl.get(50+j));
					}
					else{
						valid=false;
						break;
					}
				}
				if(valid){
					
					for (int j=1;j<=offset_right;j++){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){//&&nucl_index(nucl.get(50+j))!=-1){
							//xx[9+j]=nucl_index(nucl.get(50+j));
						}
						else{
							valid=false;
							break;
						}
					}
				}*/
				
				//System.out.println(valid);
				
		/*		
				if(valid){
					
					if(nucl.get(50).equals("A")){
						for (int i=0;i<lambda.length;i++){
							for (int k=0;k<tt1.length;k++){
								//lambda[i][k]=lambda_context_product1[i][tt1[k]][3-xx[19]][3-xx[18]][3-xx[17]][3-xx[16]][3-xx[15]]*lambda_context_product2[i][tt1[k]][3-xx[14]][3-xx[13]][3-xx[12]][3-xx[11]][3-xx[10]]
								//	*lambda_context_product3[i][tt1[k]][3-xx[9]][3-xx[8]][3-xx[7]][3-xx[6]][3-xx[5]]*lambda_context_product4[i][tt1[k]][3-xx[4]][3-xx[3]][3-xx[2]][3-xx[1]][3-xx[0]];
								lambda[i][k]=lambda_context_product1[i][tt1[k]][3-nucl_index.get(60)][3-nucl_index.get(59)][3-nucl_index.get(58)][3-nucl_index.get(57)][3-nucl_index.get(56)]*lambda_context_product2[i][tt1[k]][3-nucl_index.get(55)][3-nucl_index.get(54)][3-nucl_index.get(53)][3-nucl_index.get(52)][3-nucl_index.get(51)]
										*lambda_context_product3[i][tt1[k]][3-nucl_index.get(49)][3-nucl_index.get(48)][3-nucl_index.get(47)][3-nucl_index.get(46)][3-nucl_index.get(45)]*lambda_context_product4[i][tt1[k]][3-nucl_index.get(44)][3-nucl_index.get(43)][3-nucl_index.get(42)][3-nucl_index.get(41)][3-nucl_index.get(40)];
								
							}
						}
					}
					else if(nucl.get(50).equals("C")){
						for (int i=0;i<lambda.length;i++){
							for (int k=0;k<tt2.length;k++){							
								lambda[i][k]=lambda_context_product1[i][tt2[k]][nucl_index.get(40)][nucl_index.get(41)][nucl_index.get(42)][nucl_index.get(43)][nucl_index.get(44)]*lambda_context_product2[i][tt2[k]][nucl_index.get(45)][nucl_index.get(46)][nucl_index.get(47)][nucl_index.get(48)][nucl_index.get(49)]
										*lambda_context_product3[i][tt2[k]][nucl_index.get(51)][nucl_index.get(52)][nucl_index.get(53)][nucl_index.get(54)][nucl_index.get(55)]*lambda_context_product4[i][tt2[k]][nucl_index.get(56)][nucl_index.get(57)][nucl_index.get(58)][nucl_index.get(59)][nucl_index.get(60)];
							}
						}
					}
					else if(nucl.get(50).equals("G")){
						for (int i=0;i<lambda.length;i++){
							for (int k=0;k<tt2.length;k++){
								lambda[i][k]=lambda_context_product1[i][tt2[k]][3-nucl_index.get(60)][3-nucl_index.get(59)][3-nucl_index.get(58)][3-nucl_index.get(57)][3-nucl_index.get(56)]*lambda_context_product2[i][tt2[k]][3-nucl_index.get(55)][3-nucl_index.get(54)][3-nucl_index.get(53)][3-nucl_index.get(52)][3-nucl_index.get(51)]
									*lambda_context_product3[i][tt2[k]][3-nucl_index.get(49)][3-nucl_index.get(48)][3-nucl_index.get(47)][3-nucl_index.get(46)][3-nucl_index.get(45)]*lambda_context_product4[i][tt2[k]][3-nucl_index.get(44)][3-nucl_index.get(43)][3-nucl_index.get(42)][3-nucl_index.get(41)][3-nucl_index.get(40)];
							}
						}
					}
					else if(nucl.get(50).equals("T")){
						for (int i=0;i<lambda.length;i++){
							for (int k=0;k<tt2.length;k++){
								//lambda[i][k]=lambda_context_product1[i][tt1[k]][xx[0]][xx[1]][xx[2]][xx[3]][xx[4]]*lambda_context_product2[i][tt1[k]][xx[5]][xx[6]][xx[7]][xx[8]][xx[9]]
								//	*lambda_context_product3[i][tt1[k]][xx[10]][xx[11]][xx[12]][xx[13]][xx[14]]*lambda_context_product4[i][tt1[k]][xx[15]][xx[16]][xx[17]][xx[18]][xx[19]];
								lambda[i][k]=lambda_context_product1[i][tt1[k]][nucl_index.get(40)][nucl_index.get(41)][nucl_index.get(42)][nucl_index.get(43)][nucl_index.get(44)]*lambda_context_product2[i][tt1[k]][nucl_index.get(45)][nucl_index.get(46)][nucl_index.get(47)][nucl_index.get(48)][nucl_index.get(49)]
										*lambda_context_product3[i][tt1[k]][nucl_index.get(51)][nucl_index.get(52)][nucl_index.get(53)][nucl_index.get(54)][nucl_index.get(55)]*lambda_context_product4[i][tt1[k]][nucl_index.get(56)][nucl_index.get(57)][nucl_index.get(58)][nucl_index.get(59)][nucl_index.get(60)];
							
							}
						}
					}
					
				}
				else{
					if(nucl.get(50).equals("A")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][1];
							lambda[i][1]*=lambda_type.get(i)[0][1];
							lambda[i][2]*=lambda_type.get(i)[0][1];
							for (int j=-1;j>=-offset_left;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){
											lambda[i][k]*=lambda_context.get(i)[9-j][tt1[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10-j][tt1[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}		
							}
						}
					}
					else if(nucl.get(50).equals("C")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][0];
							lambda[i][1]*=lambda_type.get(i)[0][0];
							lambda[i][2]*=lambda_type.get(i)[0][0];
							
							for (int j=-1;j>=-offset_left;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10+j][tt2[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){
											lambda[i][k]*=lambda_context.get(i)[9+j][tt2[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
						}
					}
					else if(nucl.get(50).equals("G")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][0];
							lambda[i][1]*=lambda_type.get(i)[0][0];
							lambda[i][2]*=lambda_type.get(i)[0][0];
							
							for (int j=-1;j>=-offset_left;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){		
											lambda[i][k]*=lambda_context.get(i)[9-j][tt2[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;	
								}
							}
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt2.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10-j][tt2[k]][3-nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
								
							}
						}
					}
					else if(nucl.get(50).equals("T")){
						for (int i=0;i<lambda.length;i++){
							lambda[i][0]=lambda_type.get(i)[1][0];
							lambda[i][1]=lambda_type.get(i)[1][1];
							lambda[i][2]=lambda_type.get(i)[1][2];
							
							lambda[i][0]*=lambda_type.get(i)[0][1];
							lambda[i][1]*=lambda_type.get(i)[0][1];
							lambda[i][2]*=lambda_type.get(i)[0][1];
							
							
							for (int j=-1;j>=-offset_left;j--){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){
											lambda[i][k]*=lambda_context.get(i)[10+j][tt1[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)){
									if(nucl_index(nucl.get(50+j))!=-1){
										for (int k=0;k<tt1.length;k++){	
											lambda[i][k]*=lambda_context.get(i)[9+j][tt1[k]][nucl_index(nucl.get(50+j))];
											//x++;
										}
									}
								}
								else{
									break;
								}
							}
							
						}
					}
				}
				for (int k=0;k<lambda.length;k++){
					for (int l=0;l<lambda[k].length;l++){
						if(lambda[k][l]>1000){
							lambda[k][l]=1000;
						}
					}
				}
			*/	
				
				
				
				
			}
			catch(Exception e){
				StackTraceElement[] aa=e.getStackTrace();
				for (int i=0;i<aa.length;i++){
					System.out.println(i+"	"+aa[i].getLineNumber());
				}
				System.out.println(e);
			}
			
			
		}
		
		
		
		
		public static int nucl_index(String s){
			Integer ii=table_nucl.get(s);
			if(ii!=null){
				return ii.intValue();
			}
			else{
				return -1;
			}
		}
	
	
	public static double[][] product(double[][] x, ArrayList<Integer>[] index, double[][] y){
		double[][] z=new double[x.length][y[0].length];
		for (int i=0;i<x.length;i++){
			for (int j=0;j<y[0].length;j++){
				
				
				for (int k=0;k<index[i].size();k++){
					z[i][j]+=x[i][index[i].get(k)]*y[index[i].get(k)][j];
				}
				
//				for (int k=0;k<y.length;k++){
//					z[i][j]+=x[i][k]*y[k][j];
//				}
			}
		}
		return z;
	}
		
	public static double[][][] freq(int[][][] a){
		double[][][] b=new double[a.length][a[0].length][a[0][0].length];
		for (int i=0;i<a.length;i++){
			for (int j=0;j<a[i].length;j++){
				int sum=0;
				for (int k=0;k<a[i][j].length;k++){
					sum+=a[i][j][k];
				}
				if(sum>=50){
					for (int k=0;k<a[i][j].length;k++){
						if(a[i][j][k]==0){
							a[i][j][k]++;
							sum++;
						}
					}
					for (int k=0;k<a[i][j].length;k++){
						b[i][j][k]=(double)(a[i][j][k])/(double)(sum);
					}
				}
				else{
					
					for (int k=0;k<a[i][j].length;k++){
						b[i][j][k]=(double)((double)(50-sum)/(double)(a[i][j].length)+a[i][j][k])/(double)(50);
								
					}
					
					
					/*for (int k=0;k<a[i][j].length;k++){
						b[i][j][k]=1.0/(double)(a[i][j].length);//(double)(a[i][j][k])/(double)(sum);
					}*/
				}
			}
		}
		return b;
	}
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	public static int index_gene(String s, ArrayList<Gene> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).name.equals(s)){
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
	
	public static int min(ArrayList<int[]> a){
		int min=1000000000;
		for (int i=0;i<a.size();i++){
			if(a.get(i)[0]<min){
				min=a.get(i)[0];
			}
		}
		return min;
	}
	
	public static double gamma(double x){
		return Math.exp(Gamma.logGamma(x));
	}
	
	public static double besselk(double n, double x){
		return Bessel.k(x, n, false);
		
	}
	
	public static double log_besselk(double n, double x){
		
		
		
		double xx=Math.log(Bessel.k(x, n, false));
		//System.out.println(n+"	"+x+"	"+table_bessel[(int)(10*(n+500.0+0.051))][(int)(10*(x+0.051))]+"	"+xx);
		
		if(Double.isInfinite(xx)||Double.isNaN(xx)){
			if(0<=(int)(10*(n+500.0+0.051))&&(int)(10*(n+500.0+0.051))<table_bessel.length&&0<=(int)(10*(x+0.051))&&(int)(10*(x+0.051))<table_bessel[0].length){
				//System.out.println(n+"	"+x);
				return table_bessel[(int)(10*(n+500.0+0.051))][(int)(10*(x+0.051))];
				//return table_bessel[(int)((n+0.05)*10)][(int)((x+0.05)*10)];
			}
			else{
				return Double.NaN;
			}
			//System.out.println("NaN	"+n+"	"+x);
		}
		else{
			return xx;
		}
		
	}
	/*
	public static double log_besselk(double n, double x){
		//System.out.println(n+"	"+x+"	"+Bessel.k(x, n, false));
		
		double xx=Math.log(Bessel.k(x, n, false));
		if(Double.isFinite(xx)&&!Double.isNaN(xx)){
		//	return xx;
		}
		
		ArrayList<double[]> array=new ArrayList<double[]>();
		for (double x1=x;x1<=100;x1*=1.01){
			if(Double.isFinite(Bessel.k(x1, n, false))){
				array.add(new double[]{Math.log(x1),Math.log(Bessel.k(x1, n, false))});
			}
			if(array.size()>3){
				break;
			}
		}
		if(array.size()<2){
			//System.out.println("Bessel	"+Double.NaN);
			return Double.NaN;
		}
		
		double avg1=0;
		double avg2=0;
		for(int i=0;i<array.size();i++){
			avg1+=array.get(i)[0];
			avg2+=array.get(i)[1];
		}
		
		double sum11=0;
		double sum12=0;
		//double sum22=0;
		for(int i=0;i<array.size();i++){
			sum11+=(array.get(i)[0]-avg1)*(array.get(i)[0]-avg1);
			sum12+=(array.get(i)[0]-avg1)*(array.get(i)[1]-avg2);
		//	sum22+=(array.get(i)[1]-avg2)*(array.get(i)[1]-avg2);
		}
		double m=sum12/sum11;
		double nn=avg2-m*avg1;
		System.out.println("Bessel	"+( Math.log(x)*m+nn)+"	"+Math.log(Bessel.k(x, n, false)));
		
		return Math.log(x)*m+nn;
	}*/
	
	
	public static double pofs(double s, double L, double[] p, int modC){
		double a=0;
		double b=0;
		double t=0;
		double w=0;
		double g=0;
		double d=0;
		
		if(modC==1||modC==2){
			a=p[0];
			b=p[1];
		}
		else if(modC==3||modC==4){
			a=p[0];
			b=p[1];
			t=p[2];
			w=p[3];
		}
		else if(modC==5||modC==6){
			a=p[0];
			b=p[1];
			g=p[2];
			d=p[3];
			w=p[4];
		}
		
		
		
		if (modC==1){
			/*************** lambda ~ Gamma:*/
			
			//return (L*b)**s * (1 + L*b)**(-s-a) * math.gamma(s + a) / (math.gamma(s+1) * math.gamma(a))
			return Math.pow((L*b),s) * Math.pow((1 + L*b),(-s-a)) * gamma(s + a) / (gamma(s+1) * gamma(a));
			
			
		}
		else if (modC==2){
			/*************** lambda ~ IG:*/
			
			//return 2. * mp.exp( ((s + a)/2.)*math.log(L*b) + mp.log(mp.besselk(-s + a, 2*np.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
			return 2.0 * Math.exp( ((s + a)/2.0)*Math.log(L*b) + log_besselk(-s + a, 2*Math.sqrt(L*b)) - Gamma.logGamma(s+1) - Gamma.logGamma(a) );
			
			
		}
		else if (modC==3){
			/************** lambda ~ w * Exp + (1-w) * Gamma:*/
			
			//return np.exp( np.log(w) + s*np.log(L) + np.log(t) + (-1 - s)*np.log(L + t) ) + np.exp( np.log(1.-w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) )
			return Math.exp( Math.log(w) + s*Math.log(L) + Math.log(t) + (-1 - s)*Math.log(L + t) ) + Math.exp( Math.log(1.-w) + s*Math.log(L*b) + (-s-a)*Math.log(1 + L*b) + Gamma.logGamma(s + a) - Gamma.logGamma(s+1) - Gamma.logGamma(a) );
			
			
		}
		else if (modC==4){/// same
			/*************** lambda ~ w * Exp + (1-w) * InvGamma:*/
			
			return (w * t * Math.exp( s*Math.log(L) + (-1 - s)*Math.log(L + t) )) + Math.exp( Math.log(1.-w) + Math.log(2.) + ((s + a)/2.)*Math.log(L*b) + log_besselk(-s + a, 2*Math.sqrt(L*b)) - Gamma.logGamma(s+1) - Gamma.logGamma(a) );
			
			
		}
		else if (modC==5){
			/*************** lambda ~ w * Gamma + (1-w) * Gamma (Gamma mixture model):*/
			
			//return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + np.exp( np.log(1.-w) + s*np.log(L*d) + (-s-g)*np.log(1 + L*d) + sp.gammaln(s + g) - sp.gammaln(s+1) - sp.gammaln(g) )
			return Math.exp( Math.log(w) + s*Math.log(L*b) + (-s-a)*Math.log(1 + L*b) + Gamma.logGamma(s + a) - Gamma.logGamma(s+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1.-w) + s*Math.log(L*d) + (-s-g)*Math.log(1 + L*d) + Gamma.logGamma(s + g) - Gamma.logGamma(s+1) - Gamma.logGamma(g) );
			
			
		}
		else if (modC==6){//same
			/*************** lambda ~ w * Gamma + (1-w) * InvGamma (mixture model):*/
			
			return Math.exp( Math.log(w) + s*Math.log(L*b) + (-s-a)*Math.log(1 + L*b) + Gamma.logGamma(s + a) - Gamma.logGamma(s+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1.-w) + Math.log(2.) + ((s + g)/2.)*Math.log(L*d) + log_besselk(-s + g, 2*Math.sqrt(L*d)) - Gamma.logGamma(s+1) - Gamma.logGamma(g) );
			
		}
		else{
			return Double.NaN;
		}
	}
	
	
	public static double pofx_given_s(double x, double s, double L, double r, double[] p, int modC ){
		double a=0;
		double b=0;
		double t=0;
		double w=0;
		double g=0;
		double d=0;
		
		if(modC==1||modC==2){
			a=p[0];
			b=p[1];
		}
		else if(modC==3||modC==4){
			a=p[0];
			b=p[1];
			t=p[2];
			w=p[3];
		}
		else if(modC==5||modC==6){
			a=p[0];
			b=p[1];
			g=p[2];
			d=p[3];
			w=p[4];
		}
		
		
		
		if (modC==1){
			/*************** lambda ~ Gamma:*/
			//return np.exp( x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) / pofs(s, L);
			return Math.exp( x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==2){//same
			/*************** lambda ~ IG:*/
			
			//return mp.exp( np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (1/2. * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + mp.log(mp.besselk(s + x - a, 2*math.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) / pofs(s, L, thr);
			return Math.exp( Math.log(2) + (s + x)*Math.log(L) + x*Math.log(r) + (1/2. * (-s - x + a))*Math.log((L*(1 + r))/b) + a*Math.log(b) + log_besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b)) -  Gamma.logGamma(s+1) -  Gamma.logGamma(x+1) -  Gamma.logGamma(a) ) / pofs(s, L , p , modC);
				
			
		}
		else if (modC==3){
			/************** lambda ~ w * Exp + (1-w) * Gamma:*/
			//return ( np.exp( np.log(w) + (s + x)*np.log(L) + x*np.log(r) + np.log(t) + (-1 - s - x)*np.log(L + L*r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s+1) - sp.gammaln(x+1) ) + np.exp( np.log(1-w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) ) / pofs(s, L);
			return ( Math.exp( Math.log(w) + (s + x)*Math.log(L) + x*Math.log(r) + Math.log(t) + (-1 - s - x)*Math.log(L + L*r + t) + Gamma.logGamma(1 + s + x) -Gamma.logGamma(s+1) - Gamma.logGamma(x+1) ) + Math.exp( Math.log(1-w) + x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==4){//same
			/*************** lambda ~ w * Exp + (1-w) * InvGamma:*/
			
			//System.out.println((s + x - a)+"	"+(2*Math.sqrt(L*(1 + r)*b))+"	"+besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b)));
			//System.out.println(Math.log(w) +"	"+ (s + x)*Math.log(L) +"	"+  x*Math.log(r) +"	"+  Math.log(t) +"	"+  (-1 - s - x)*Math.log(L + L*r + t) +"	"+  Gamma.logGamma(1 + s + x) +"	"+  Gamma.logGamma(s+1) +"	"+  Gamma.logGamma(x+1)+"	X	"+Math.log(1.-w) +"	"+  Math.log(2) +"	"+  (s + x)*Math.log(L) +"	"+  x*Math.log(r) +"	"+  (0.5 * (-s - x + a))*Math.log((L*(1 + r))/b) +"	"+  a*Math.log(b) +"	"+  Math.log(besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b))) +"	"+  Gamma.logGamma(s+1) +"	"+  Gamma.logGamma(x+1) +"	"+  Gamma.logGamma(a));
			//return ( np.exp( np.log(w) + (s + x)*np.log(L) + x*np.log(r) + np.log(t) + (-1 - s - x)*np.log(L + L*r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s+1) - sp.gammaln(x+1) ) + mp.exp( np.log(1.-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + mp.log(mp.besselk(s + x - a, 2*np.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) ) / pofs(s, L, thr);
			return ( Math.exp( Math.log(w) + (s + x)*Math.log(L) + x*Math.log(r) + Math.log(t) + (-1 - s - x)*Math.log(L + L*r + t) + Gamma.logGamma(1 + s + x) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) ) + Math.exp( Math.log(1.-w) + Math.log(2) + (s + x)*Math.log(L) + x*Math.log(r) + (0.5 * (-s - x + a))*Math.log((L*(1 + r))/b) + a*Math.log(b) + log_besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b)) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==5){
			/*************** lambda ~ w * Gamma + (1-w) * Gamma (Gamma mixture model):*/
			
			//return ( np.exp( np.log(w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) + np.exp( np.log(1-w) + x*np.log(r) + (s + x)*np.log(L*d) + (-s - x - g)*np.log(1 + L*(1 + r)*d) + sp.gammaln(s + x + g) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(g) ) ) / pofs(s, L);
			return ( Math.exp( Math.log(w) + x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + x*Math.log(r) + (s + x)*Math.log(L*d) + (-s - x - g)*Math.log(1 + L*(1 + r)*d) + Gamma.logGamma(s + x + g) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(g) ) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==6){//same
			/*************** lambda ~ w * Gamma + (1-w) * InvGamma (mixture model):*/
			
			//return ( np.exp( np.log(w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) + mp.exp( np.log(1-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + g))*np.log((L*(1 + r))/d) + g*np.log(d) + mp.log(mp.besselk(s + x - g, 2*np.sqrt(L*(1 + r)*d))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(g) ) ) / pofs(s, L, thr);
			return ( Math.exp( Math.log(w) + x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + Math.log(2) + (s + x)*Math.log(L) + x*Math.log(r) + (0.5 * (-s - x + g))*Math.log((L*(1 + r))/d) + g*Math.log(d) + log_besselk(s + x - g, 2*Math.sqrt(L*(1 + r)*d)) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(g) ) ) / pofs(s, L, p , modC);
			
			
		}
		else{
			return Double.NaN;
		}
	}
	
	
	public static double pofsX(double s, double L, double[] p, int modC){
		double a=0;
		double b=0;
		double t=0;
		double w=0;
		double g=0;
		double d=0;
		
		if(modC==1||modC==2){
			a=p[0];
			b=p[1];
		}
		else if(modC==3||modC==4){
			a=p[0];
			b=p[1];
			t=p[2];
			w=p[3];
		}
		else if(modC==5||modC==6){
			a=p[0];
			b=p[1];
			g=p[2];
			d=p[3];
			w=p[4];
		}
		
		
		
		if (modC==1){
			/*************** lambda ~ Gamma:*/
			
			//return (L*b)**s * (1 + L*b)**(-s-a) * math.gamma(s + a) / (math.gamma(s+1) * math.gamma(a))
			return Math.pow((L*b),s) * Math.pow((1 + L*b),(-s-a)) * gamma(s + a) / (gamma(s+1) * gamma(a));
			
			
		}
		else if (modC==2){
			/*************** lambda ~ IG:*/
			
			//return 2. * mp.exp( ((s + a)/2.)*math.log(L*b) + mp.log(mp.besselk(-s + a, 2*np.sqrt(L*b))) - sp.gammaln(s+1) - sp.gammaln(a) )
			return 2.0 * Math.exp( ((s + a)/2.0)*Math.log(L*b) + Math.log(besselk(-s + a, 2*Math.sqrt(L*b))) - Gamma.logGamma(s+1) - Gamma.logGamma(a) );
			
			
		}
		else if (modC==3){
			/************** lambda ~ w * Exp + (1-w) * Gamma:*/
			
			//return np.exp( np.log(w) + s*np.log(L) + np.log(t) + (-1 - s)*np.log(L + t) ) + np.exp( np.log(1.-w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) )
			return Math.exp( Math.log(w) + s*Math.log(L) + Math.log(t) + (-1 - s)*Math.log(L + t) ) + Math.exp( Math.log(1.-w) + s*Math.log(L*b) + (-s-a)*Math.log(1 + L*b) + Gamma.logGamma(s + a) - Gamma.logGamma(s+1) - Gamma.logGamma(a) );
			
			
		}
		else if (modC==4){/// same
			/*************** lambda ~ w * Exp + (1-w) * InvGamma:*/
			
			return (w * t * Math.exp( s*Math.log(L) + (-1 - s)*Math.log(L + t) )) + Math.exp( Math.log(1.-w) + Math.log(2.) + ((s + a)/2.)*Math.log(L*b) + Math.log(besselk(-s + a, 2*Math.sqrt(L*b))) - Gamma.logGamma(s+1) - Gamma.logGamma(a) );
			
			
		}
		else if (modC==5){
			/*************** lambda ~ w * Gamma + (1-w) * Gamma (Gamma mixture model):*/
			
			//return np.exp( np.log(w) + s*np.log(L*b) + (-s-a)*np.log(1 + L*b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a) ) + np.exp( np.log(1.-w) + s*np.log(L*d) + (-s-g)*np.log(1 + L*d) + sp.gammaln(s + g) - sp.gammaln(s+1) - sp.gammaln(g) )
			return Math.exp( Math.log(w) + s*Math.log(L*b) + (-s-a)*Math.log(1 + L*b) + Gamma.logGamma(s + a) - Gamma.logGamma(s+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1.-w) + s*Math.log(L*d) + (-s-g)*Math.log(1 + L*d) + Gamma.logGamma(s + g) - Gamma.logGamma(s+1) - Gamma.logGamma(g) );
			
			
		}
		else if (modC==6){//same
			/*************** lambda ~ w * Gamma + (1-w) * InvGamma (mixture model):*/
			
			return Math.exp( Math.log(w) + s*Math.log(L*b) + (-s-a)*Math.log(1 + L*b) + Gamma.logGamma(s + a) - Gamma.logGamma(s+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1.-w) + Math.log(2.) + ((s + g)/2.)*Math.log(L*d) + Math.log(besselk(-s + g, 2*Math.sqrt(L*d))) - Gamma.logGamma(s+1) - Gamma.logGamma(g) );
			
		}
		else{
			return Double.NaN;
		}
	}
	
	
	public static double pofx_given_sX(double x, double s, double L, double r, double[] p, int modC ){
		double a=0;
		double b=0;
		double t=0;
		double w=0;
		double g=0;
		double d=0;
		
		if(modC==1||modC==2){
			a=p[0];
			b=p[1];
		}
		else if(modC==3||modC==4){
			a=p[0];
			b=p[1];
			t=p[2];
			w=p[3];
		}
		else if(modC==5||modC==6){
			a=p[0];
			b=p[1];
			g=p[2];
			d=p[3];
			w=p[4];
		}
		
		
		
		if (modC==1){
			/*************** lambda ~ Gamma:*/
			//return np.exp( x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) / pofs(s, L);
			return Math.exp( x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==2){//same
			/*************** lambda ~ IG:*/
			
			//return mp.exp( np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (1/2. * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + mp.log(mp.besselk(s + x - a, 2*math.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) / pofs(s, L, thr);
			return Math.exp( Math.log(2) + (s + x)*Math.log(L) + x*Math.log(r) + (1/2. * (-s - x + a))*Math.log((L*(1 + r))/b) + a*Math.log(b) + Math.log(besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b))) -  Gamma.logGamma(s+1) -  Gamma.logGamma(x+1) -  Gamma.logGamma(a) ) / pofs(s, L , p , modC);
				
			
		}
		else if (modC==3){
			/************** lambda ~ w * Exp + (1-w) * Gamma:*/
			//return ( np.exp( np.log(w) + (s + x)*np.log(L) + x*np.log(r) + np.log(t) + (-1 - s - x)*np.log(L + L*r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s+1) - sp.gammaln(x+1) ) + np.exp( np.log(1-w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) ) / pofs(s, L);
			return ( Math.exp( Math.log(w) + (s + x)*Math.log(L) + x*Math.log(r) + Math.log(t) + (-1 - s - x)*Math.log(L + L*r + t) + Gamma.logGamma(1 + s + x) -Gamma.logGamma(s+1) - Gamma.logGamma(x+1) ) + Math.exp( Math.log(1-w) + x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==4){//same
			/*************** lambda ~ w * Exp + (1-w) * InvGamma:*/
			
			//System.out.println((s + x - a)+"	"+(2*Math.sqrt(L*(1 + r)*b))+"	"+besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b)));
			//System.out.println(Math.log(w) +"	"+ (s + x)*Math.log(L) +"	"+  x*Math.log(r) +"	"+  Math.log(t) +"	"+  (-1 - s - x)*Math.log(L + L*r + t) +"	"+  Gamma.logGamma(1 + s + x) +"	"+  Gamma.logGamma(s+1) +"	"+  Gamma.logGamma(x+1)+"	X	"+Math.log(1.-w) +"	"+  Math.log(2) +"	"+  (s + x)*Math.log(L) +"	"+  x*Math.log(r) +"	"+  (0.5 * (-s - x + a))*Math.log((L*(1 + r))/b) +"	"+  a*Math.log(b) +"	"+  Math.log(besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b))) +"	"+  Gamma.logGamma(s+1) +"	"+  Gamma.logGamma(x+1) +"	"+  Gamma.logGamma(a));
			//return ( np.exp( np.log(w) + (s + x)*np.log(L) + x*np.log(r) + np.log(t) + (-1 - s - x)*np.log(L + L*r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s+1) - sp.gammaln(x+1) ) + mp.exp( np.log(1.-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + a))*np.log((L*(1 + r))/b) + a*np.log(b) + mp.log(mp.besselk(s + x - a, 2*np.sqrt(L*(1 + r)*b))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) ) / pofs(s, L, thr);
			return ( Math.exp( Math.log(w) + (s + x)*Math.log(L) + x*Math.log(r) + Math.log(t) + (-1 - s - x)*Math.log(L + L*r + t) + Gamma.logGamma(1 + s + x) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) ) + Math.exp( Math.log(1.-w) + Math.log(2) + (s + x)*Math.log(L) + x*Math.log(r) + (0.5 * (-s - x + a))*Math.log((L*(1 + r))/b) + a*Math.log(b) + Math.log(besselk(s + x - a, 2*Math.sqrt(L*(1 + r)*b))) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==5){
			/*************** lambda ~ w * Gamma + (1-w) * Gamma (Gamma mixture model):*/
			
			//return ( np.exp( np.log(w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) + np.exp( np.log(1-w) + x*np.log(r) + (s + x)*np.log(L*d) + (-s - x - g)*np.log(1 + L*(1 + r)*d) + sp.gammaln(s + x + g) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(g) ) ) / pofs(s, L);
			return ( Math.exp( Math.log(w) + x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + x*Math.log(r) + (s + x)*Math.log(L*d) + (-s - x - g)*Math.log(1 + L*(1 + r)*d) + Gamma.logGamma(s + x + g) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(g) ) ) / pofs(s, L, p , modC);
			
		}
		else if (modC==6){//same
			/*************** lambda ~ w * Gamma + (1-w) * InvGamma (mixture model):*/
			
			//return ( np.exp( np.log(w) + x*np.log(r) + (s + x)*np.log(L*b) + (-s - x - a)*np.log(1 + L*(1 + r)*b) + sp.gammaln(s + x + a) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(a) ) + mp.exp( np.log(1-w) + np.log(2) + (s + x)*np.log(L) + x*np.log(r) + (0.5 * (-s - x + g))*np.log((L*(1 + r))/d) + g*np.log(d) + mp.log(mp.besselk(s + x - g, 2*np.sqrt(L*(1 + r)*d))) - sp.gammaln(s+1) - sp.gammaln(x+1) - sp.gammaln(g) ) ) / pofs(s, L, thr);
			return ( Math.exp( Math.log(w) + x*Math.log(r) + (s + x)*Math.log(L*b) + (-s - x - a)*Math.log(1 + L*(1 + r)*b) + Gamma.logGamma(s + x + a) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(a) ) + Math.exp( Math.log(1-w) + Math.log(2) + (s + x)*Math.log(L) + x*Math.log(r) + (0.5 * (-s - x + g))*Math.log((L*(1 + r))/d) + g*Math.log(d) + Math.log(besselk(s + x - g, 2*Math.sqrt(L*(1 + r)*d))) - Gamma.logGamma(s+1) - Gamma.logGamma(x+1) - Gamma.logGamma(g) ) ) / pofs(s, L, p , modC);
			
			
		}
		else{
			return Double.NaN;
		}
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
	
	private static class Gene{
		String name="";
		ArrayList<int[]> coord=new ArrayList<int[]>();
		int start=-1;
		int end=-1;
		
		double cov=0;
		double cov_syn=0;
		int[] count=null;
		int[] count_syn=null;
		int[] count_all=null;
		int[] count_destruct=null;
		
		double[] ell_s=null;
		double[] ell_x=null;
		double[] sobs=null;
		double[] xobs=null;
		
		double[] sign_combined=null;
		double[] sign_combined_uniform=null;
		double[] sign_destruct=null;
		double[] sign_hotspot=null;
		double[] sign_complete=null;
		double[] sign_complete_uniform=null;
		double[] sign_vector_syn=null;
		double[] sign_hotspot_syn=null;
		double[] sign_complete_syn=null;
		
		public Gene(String name, int start, int end, int no_entities){
			this.name=name;
			coord.add(new int[]{start,end});
			
			this.count=new int[no_entities];
			this.count_syn=new int[no_entities];
			this.count_all=new int[no_entities];
			this.count_destruct=new int[no_entities];
			
			this.ell_s=new double[no_entities];
			this.ell_x=new double[no_entities];
			this.sobs=new double[no_entities];
			this.xobs=new double[no_entities];
			
			this.sign_combined=new double[no_entities];
			this.sign_combined_uniform=new double[no_entities];
			this.sign_destruct=new double[no_entities];
			this.sign_hotspot=new double[no_entities];
			this.sign_complete=new double[no_entities];
			this.sign_complete_uniform=new double[no_entities];
			this.sign_vector_syn=new double[no_entities];
			this.sign_hotspot_syn=new double[no_entities];
			this.sign_complete_syn=new double[no_entities];
			
			for (int i=0;i<no_entities;i++){
				this.sign_combined[i]=1;
				this.sign_combined_uniform[i]=1;
				this.sign_destruct[i]=1;
				this.sign_hotspot[i]=1;
				this.sign_complete[i]=1;
				this.sign_complete_uniform[i]=1;
				this.sign_vector_syn[i]=1;
				this.sign_hotspot_syn[i]=1;
				this.sign_complete_syn[i]=1;
				
				this.ell_s[i]=Double.NaN;
				this.sobs[i]=Double.NaN;
				this.ell_x[i]=Double.NaN;
				this.xobs[i]=Double.NaN;
				
			}
			
		}
		public boolean contains (int pos){
			for (int i=0;i<coord.size();i++){
				if(coord.get(i)[0]<=pos&&pos<=coord.get(i)[1]){
					return true;
				}
			}
			return false;
		}
		
		public GeneSmall small (int k){
			GeneSmall g=new GeneSmall();
			
			g.name=this.name;
			g.coord=this.coord;
			g.cov=this.cov;
			g.cov_syn=this.cov_syn;
			g.count=this.count[k];
			g.count_syn=this.count_syn[k];
			g.count_all=this.count_all[k];
			g.count_destruct=this.count_destruct[k];
			
			g.sign_combined=this.sign_combined[k];
			g.sign_destruct=this.sign_destruct[k];
			g.sign_hotspot=this.sign_hotspot[k];
			g.sign_complete=this.sign_complete[k];
			g.sign_vector_syn=this.sign_vector_syn[k];
			g.sign_hotspot_syn=this.sign_hotspot_syn[k];
			g.sign_complete_syn=this.sign_complete_syn[k];
			
			return g;
		}
		
		public GeneSmall small_uniform (int k){
			GeneSmall g=new GeneSmall();
			
			g.name=this.name;
			g.coord=this.coord;
			g.cov=this.cov;
			g.cov_syn=this.cov_syn;
			g.count=this.count[k];
			g.count_syn=this.count_syn[k];
			g.count_all=this.count_all[k];
			g.count_destruct=this.count_destruct[k];
			
			g.sign_combined=this.sign_combined[k];
			g.sign_combined=this.sign_combined_uniform[k];
			g.sign_destruct=this.sign_destruct[k];
			g.sign_hotspot=this.sign_hotspot[k];
			g.sign_complete=this.sign_complete[k];
			g.sign_complete=this.sign_complete_uniform[k];
			g.sign_vector_syn=this.sign_vector_syn[k];
			g.sign_hotspot_syn=this.sign_hotspot_syn[k];
			g.sign_complete_syn=this.sign_complete_syn[k];
			
			return g;
		}
	}
	
	private static class GeneSmall{
		String name="";
		ArrayList<int[]> coord=new ArrayList<int[]>();
		int start=-1;
		int end=-1;
		
		double cov=0;
		double cov_syn=0;
		
		int count=0;
		int count_syn=0;
		int count_all=0;
		int count_destruct=0;
		
		double sign_combined=1;
		//double sign_combined_uniform=1;
		double sign_destruct=1;
		double sign_hotspot=1;
		double sign_complete=1;
		//double sign_complete_uniform=1;
		double sign_vector_syn=1;
		double sign_hotspot_syn=1;
		double sign_complete_syn=1;
		double fdr=1;
		double fdr_syn=1;
		
		
		public boolean contains (int pos){
			for (int i=0;i<coord.size();i++){
				if(coord.get(i)[0]<=pos&&pos<=coord.get(i)[1]){
					return true;
				}
			}
			return false;
		}
	}
}
