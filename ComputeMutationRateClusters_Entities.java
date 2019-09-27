/************************************************************           
 * MutPanning 												*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2019 									*   
 *															*   
 * License:		BSD-3-Clause open source license			*  
 *															*   
 * Summary: This script computes the mutation rate of each	*
 * 			cluster. These clusters were dervied in the 	*
 * 			previous steps and contain samples with 		*
 * 			similar passenger mutation distribution patterns*
 * 			These mutation rates are then summed up to		*
 * 			compare the syn. vs. nonsyn. target size of		*
 * 			each gene.										*
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




public class ComputeMutationRateClusters_Entities {
	static ArrayList<double[][][]> lambda_context=new ArrayList<double[][][]>();
	static ArrayList<double[][]> lambda_type=new ArrayList<double[][]>();
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
	
	//static double[][][][][][][][] lambda_context_product5=null;
	
	
	static ArrayList<Integer> pos=null;
	static ArrayList<String> nucl=null;
	static ArrayList<Integer> nucl_index=null;
	static ArrayList<Double> coverage=null;
	static ArrayList<int[]> label=null;
	static ArrayList<int[][]> count=null;
	
	
	
	static String[] entities=null;
	static String[] index_header_samples={"ID","Sample","Cohort"};
	static int no_clusters=-1;
	static ArrayList<Gene>[] genes=new ArrayList[chr.length];
	//static ArrayList<Integer>[] list=new ArrayList[0];
	static Hashtable<Integer,Integer> table_entity=null;
	
	//static String file_out="";
	static String file_annotation="";
	static String file_align="";
	static String file_signatures="";
	static String file_reference="";
	static String file_clusters="";
	static String file_type="";
	static String file_samples="";
	static String file_genes="";
	static String file_out2="";
	
	
	static double[][] weights=null;
	
	/*
	 * argument0: root file
	 * 
	 */
	
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
	
	static double[][][] lambda=null;
	//static int ii=-1;
	public static void main(String[] args){
		
		//file_out=args[0]+"MutationRateClusters/Lambda_Chr";
		file_annotation=args[2]+"AnnotationHg19/Annotation_chr";
		file_align=args[0]+"AlignHg19/AlignHg19Chr";
		file_signatures=args[0]+"ClusteringComplete/ClusteringComplete_Affinity.txt";
		file_reference=args[2]+"FileReferenceCount.txt";
		file_clusters=args[0]+"ClusteringComplete/ClusteringComplete_Samples.txt";
		file_type=args[0]+"AffinityCounts/TypeCount.txt";
		file_samples=args[1];
		file_genes=args[2]+"Exons_Hg19.txt";
		file_out2=args[0]+"CBASE/CountsRaw/Count";
		
		
		if (!new File(args[0]+"CBASE/CountsRaw/").exists()){
			new File(args[0]+"CBASE/CountsRaw/").mkdirs();
		}
		
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
			
			//each cluster delivers a "signature" for the distribution of passenger
			//mutations in that cluster. read that signature
			ArrayList<double[][][]> frequency=new ArrayList<double[][][]>();
			FileInputStream in=new FileInputStream(file_signatures);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
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
			
			double[][][] lambda_context_product5=new double[lambda_context.size()][6][4096];

			int[] ttt={0,1,2,0,2,1};
			
			for (int a=0;a<lambda_context.size();a++){
				
				for (int k=0;k<6;k++){
					int[] x=new int[6];
					for (x[0]=0;x[0]<4;x[0]++){
						for (x[1]=0;x[1]<4;x[1]++){
							for (x[2]=0;x[2]<4;x[2]++){
								for (x[3]=0;x[3]<4;x[3]++){
									for (x[4]=0;x[4]<4;x[4]++){
										for (x[5]=0;x[5]<4;x[5]++){
											
											int n=1;
											int index=0;
											for (int l=0;l<x.length;l++){
												index+=x[l]*n;
												n*=4;
											}
											
											lambda_context_product5[a][k][index]=lambda_type.get(a)[0][k/3]*lambda_type.get(a)[1][ttt[k]];
											
											for (int l=0;l<x.length;l++){
												lambda_context_product5[a][k][index]*=lambda_context.get(a)[l+7][k][x[l]];
											}
									
										
										}
										
									}
								}
							}
						}
					}
				}
				
			}
			
			
			
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
			}
			
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
			
			lambda=new double[weights.length][6][4096];
			for (int i=0;i<weights.length;i++){
				for (int j=0;j<lambda[i].length;j++){
					for (int k=0;k<lambda[i][j].length;k++){
						for(int l=0;l<weights[i].length;l++){
							lambda[i][j][k]+=weights[i][l]*lambda_context_product5[l][j][k];
						}
					}
				}
			}
			//ii=index("Skin",entities);
			
			
			
			
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
					genes[Integer.parseInt(t[1])-1].add(new Gene(t[0],Integer.parseInt(t[2]),Integer.parseInt(t[3])));
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
			
			table_entity=new Hashtable<Integer,Integer>();
			//list=new ArrayList[entities.length];
			//for (int i=0;i<list.length;i++){
			//	list[i]=new ArrayList<Integer>();
			//}
			
			in=new FileInputStream(file_samples);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			index_header=index_header(input.readLine().split("	"),index_header_samples);
			s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				int ii=index(t[index_header[2]],entities);
				table_entity.put(Integer.parseInt(t[index_header[0]]),ii);
				//list[ii].add(Integer.parseInt(t[index_header[0]]));
			}
			input.close();
			
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
			
			System.out.println("START");
			for (int i=0;i<chr.length;i++){
				//run(index("Y",chr),lambda);
				run(i,lambda);
			}
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
	
	//each thread computes the clusterwise local mutation rate for 1 chromosome
	//(reading the reference sequence is usually the rate limiting step)
	//private static class Subthread extends Thread{
	
	static ArrayList<String> nucl2=new ArrayList<String>();
	static ArrayList<Integer> pos2=new ArrayList<Integer>();
	static ArrayList<Integer> index2=new ArrayList<Integer>();
	static ArrayList<Double> coverage2=new ArrayList<Double>();
	static ArrayList<int[]> label2=new ArrayList<int[]>();
	static ArrayList<int[][]> count2=new ArrayList<int[][]>();
	
		
	static int[] tt1={0,1,2};
	static int[] tt2={3,5,4};
	
		//BufferedWriter output=null;
		public static void run(int c, double[][][] lambda){
			
			
			
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
			
			
			
			//product(weights,lambda)
			
			pos=new ArrayList<Integer>();
			nucl=new ArrayList<String>();
			nucl_index=new ArrayList<Integer>();
			coverage=new ArrayList<Double>();
			label=new ArrayList<int[]>();
			count=new ArrayList<int[][]>();
			
			nucl2=new ArrayList<String>();
			pos2=new ArrayList<Integer>();
			index2=new ArrayList<Integer>();
			coverage2=new ArrayList<Double>();
			label2=new ArrayList<int[]>();
			count2=new ArrayList<int[][]>();
			
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
				
				FileWriter[] out2=new FileWriter[entities.length+1];
				BufferedWriter[] output2= new BufferedWriter[entities.length+1];
				
				for (int i=0;i<entities.length;i++){
					out2[i]=new FileWriter(file_out2+entities[i]+"_Chr"+chr[c]+".txt");
					output2[i]=new BufferedWriter(out2[i]);
				}
				out2[out2.length-1]=new FileWriter(file_out2+"PanCancer"+"_Chr"+chr[c]+".txt"); 
				output2[out2.length-1]=new BufferedWriter(out2[out2.length-1]);
				
				
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
							index2.remove(i);
							pos2.remove(i);
							coverage2.remove(i);
							count2.remove(i);
							label2.remove(i);
						}
					}
					
					
					while((s=input.readLine())!=null){
						nnn++;
						//if(nnn%10000==0){
						//	System.out.println(nnn+"/"+n);
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
							//for (int i=0;i<list.length;i++){
							//	count_local[i][k]=overlap(indices[k],list[i]);
							//}
							
						}
						
						
//						if(Integer.parseInt(t[0])==19748006){
//							System.out.println(label_local[0]+"	"+label_local[1]+"	"+label_local[2]);
//							System.out.println(pos.size()+","+nucl.size()+","+coverage.size()+","+label.size()+","+count.size());
//							System.exit(0);
//						}
						
						pos.add(Integer.parseInt(t[0]));
						nucl.add(t[1]);
						nucl_index.add(nucl_index(t[1]));
						coverage.add(Double.parseDouble(t[2]));
						label.add(label_local);
						count.add(count_local);
						
						
						//when the queue is large enough, compute the mutation rate in the cnenter
						//(update function) and delete the first element
						if(pos.size()>100){
							//update_10_10();
							update_3_3();
							
							if(index2.size()==45000){
							//	return;
							}
							pos.remove(0);
							nucl.remove(0);
							nucl_index.remove(0);
							coverage.remove(0);
							label.remove(0);
							count.remove(0);
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
							//update_10_10();
							update_3_3();
							pos.remove(0);
							nucl.remove(0);
							nucl_index.remove(0);
							coverage.remove(0);
							label.remove(0);
							count.remove(0);
						}
					}
					
					double cov_gene=0;
					double cov_gene_syn=0;
					int[] count_gene=new int[entities.length+1];
					int[] count_gene_syn=new int[entities.length+1];
					double[] lambda_gene=new double[entities.length+1];
					double[] lambda_gene_syn=new double[entities.length+1];
					
					
					
					//sum up the mutation rates and mutation counts over the gene to
					//compute its target size and total no. mutations
					//note that the mutation rate and target size is determined separately
					//for synonymous and nonsynonymous mutations
					for (int i=0;i<pos2.size();i++){
						
						if(genes[c].get(nn).contains(pos2.get(i))){
							for (int j=0;j<label2.get(i).length;j++){
								if(label2.get(i)[j]!=0){
									cov_gene_syn+=coverage2.get(i);
									for (int k=0;k<count2.get(i).length;k++){
										count_gene[k]+=count2.get(i)[k][j];
										if(index2.get(i)!=-1){
											if(nucl2.get(i).equals("A")||nucl2.get(i).equals("T")){
												lambda_gene[k]+=coverage2.get(i)*lambda[k][tt2[j]][index2.get(i)];//lambda2.get(i)[k][j]; 
											}
											else if(nucl2.get(i).equals("C")||nucl2.get(i).equals("G")){
												lambda_gene[k]+=coverage2.get(i)*lambda[k][tt1[j]][index2.get(i)];//lambda2.get(i)[k][j]; 
											}
											
										}
									}
								}
								else{
									cov_gene+=coverage2.get(i);
									for (int k=0;k<count2.get(i).length;k++){
										count_gene_syn[k]+=count2.get(i)[k][j];
										if(index2.get(i)!=-1){
											if(nucl2.get(i).equals("A")||nucl2.get(i).equals("T")){
												lambda_gene_syn[k]+=coverage2.get(i)*lambda[k][tt2[j]][index2.get(i)];//lambda2.get(i)[k][j]; 
											}
											else if(nucl2.get(i).equals("C")||nucl2.get(i).equals("G")){
												lambda_gene_syn[k]+=coverage2.get(i)*lambda[k][tt1[j]][index2.get(i)];//lambda2.get(i)[k][j]; 
											}
										}
										//lambda_gene_syn[k]+=coverage2.get(i)*lambda2.get(i)[k][j];
									}
								}
							}
						}
					}
					
					//output of the gene count and target size
					//genes which have nearly no coverage are masked as they are producing artifacts
					if(cov_gene+cov_gene_syn>30){
						//for (int i=0;i<count_gene.length;i++){
//						for (int i=0;i<entities.length;i++){
//							if(entities[i].equals("Skin")){
//								System.out.println(genes[c].get(nn).name+"	"+cov_gene+"	"+cov_gene_syn+"	"+lambda_gene[i]+"	"+lambda_gene_syn[i]+"	"+count_gene[i]+"	"+count_gene_syn[i]);
//							}
//						}
						for (int i=0;i<output2.length;i++){
							output2[i].write(genes[c].get(nn).name+"	"+cov_gene+"	"+cov_gene_syn+"	"+lambda_gene[i]+"	"+lambda_gene_syn[i]+"	"+count_gene[i]+"	"+count_gene_syn[i]);
							output2[i].newLine();
						}
					}
					
					//System.out.println(nn+"	"+genes[c].get(nn).name);
				}
					//compute the last mutation rates at the end of the chr until the queue empty 
				
				input.close();
				input2.close();
				for (int i=0;i<output2.length;i++){
					output2[i].close();
				}
				//output.close();
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
		
	
		
		//the "heart" of this thread which computes the local mutation rate in the center of a queu and appends it to output
		//in brief walks from -10 to +10 around the central position (skipping 0) and forming the product over likelihoods
		
		
	
	/*	
		
		static int[] tt1={3,5,4};
		static int[] tt2={0,1,2};
		static int offset_left=10;
		static int offset_right=10;
		public static void update_10_10(){
			try{
				//int x=0;
				double[][] lambda=new double[lambda_context.size()][3];
				//int[] xx=new int[20];
				
				
				boolean valid=true;
				
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
				}
				
				//System.out.println(valid);
				
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
				
				lambda2.add(product(weights,lambda));
				pos2.add(pos.get(50));
				nucl2.add(nucl.get(50));
				coverage2.add(coverage.get(50));
				label2.add(label.get(50));
				count2.add(count.get(50));
				
				
			}
			catch(Exception e){
				StackTraceElement[] aa=e.getStackTrace();
				for (int i=0;i<aa.length;i++){
					System.out.println(i+"	"+aa[i].getLineNumber());
				}
				System.out.println(e);
			}
			
			
		}
	*/	
		
		static int[] tt1_2={3,5,4};
		static int[] tt2_2={0,1,2};
		static int offset_left_2=3;
		static int offset_right_2=3;
		public static void update_3_3(){
			try{
				boolean valid=true;
				
				int sum=0;
				if(nucl.get(50).equals("C")||nucl.get(50).equals("T")){
					int n=1;
					for (int j=-offset_left_2;j<=-1;j++){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){//nucl_index(nucl.get(50+j))!=-1){
							sum+=nucl_index(nucl.get(50+j))*n;
						}
						else{
							valid=false;
							break;
						}
						n*=4;
					}
					
					
					for (int j=1;j<=offset_right_2;j++){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){//&&nucl_index(nucl.get(50+j))!=-1){
							sum+=nucl_index(nucl.get(50+j))*n;
						}
						else{
							valid=false;
							break;
						}
						n*=4;
					}
					
				}
				else if(nucl.get(50).equals("A")||nucl.get(50).equals("G")){
					
					int n=1;
					for (int j=offset_right_2;j>=1;j--){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){//&&nucl_index(nucl.get(50+j))!=-1){
							sum+=(3-nucl_index(nucl.get(50+j)))*n;
						}
						else{
							valid=false;
							break;
						}
						n*=4;
					}
					
					for (int j=-1;j>=-offset_left_2;j--){
						if(pos.get(50)+j==pos.get(50+j)&&nucl_index.get(50+j)!=-1){//nucl_index(nucl.get(50+j))!=-1){
							sum+=(3-nucl_index(nucl.get(50+j)))*n;
						}
						else{
							valid=false;
							break;
						}
						n*=4;
					}
				}
				
				if (valid){
					index2.add(sum);
					/*
					System.out.print(pos.get(50)+"	"+nucl.get(50)+"	"+coverage.get(50));
					for (int j=0;j<3;j++){
						if(nucl.get(50).equals("A")||nucl.get(50).equals("T")){
							System.out.print("	"+lambda[ii][tt2[j]][sum]);
						}
						else if(nucl.get(50).equals("C")||nucl.get(50).equals("G")){
							System.out.print("	"+lambda[ii][tt1[j]][sum]);
						}
					}
					System.out.println();
					*/
					
					
				}
				else{
					index2.add(-1);
				}
				pos2.add(pos.get(50));
				nucl2.add(nucl.get(50));
				coverage2.add(coverage.get(50));
				label2.add(label.get(50));
				count2.add(count.get(50));
				
			}
			catch(Exception e){
				StackTraceElement[] aa=e.getStackTrace();
				for (int i=0;i<aa.length;i++){
					System.out.println(i+"	"+aa[i].getLineNumber());
				}
				System.out.println(e);
			}
			
			
		}
		
		
		
		
		
		public static void updateXX(){
			int[] tt1={3,5,4};
			int[] tt2={0,1,2};
			int offset_left=10;
			int offset_right=10;
					
			
			
			try{
				//int x=0;
				double[][] lambda=new double[lambda_context.size()][3];
				if(nucl.get(50).equals("A")){
					for (int i=0;i<lambda.length;i++){
						
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][1];
						lambda[i][1]*=lambda_type.get(i)[0][1];
						lambda[i][2]*=lambda_type.get(i)[0][1];
						//int[] tt={3,5,4};
						
						
						//for (int j=-offset_left;j<0;j++){
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
						
						/*
						int[] xx=new int[20];
						boolean valid=true;
						for (int j=-1;j>=-offset_left;j--){
							if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
								xx[10+j]=nucl_index(nucl.get(50+j));
							}
							else{
								valid=false;
								break;
							}
						}
						
						if(valid){
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									xx[9+j]=nucl_index(nucl.get(50+j));
								}
								else{
									valid=false;
									break;
								}
							}
						}
						
						if (valid){
							//double[] p=new double[3];
							for (int k=0;k<tt1.length;k++){
								
								lambda[i][k]=lambda_context_product1[i][tt1[k]][3-xx[19]][3-xx[18]][3-xx[17]][3-xx[16]][3-xx[15]]*lambda_context_product2[i][tt1[k]][3-xx[14]][3-xx[13]][3-xx[12]][3-xx[11]][3-xx[10]]
										*lambda_context_product3[i][tt1[k]][3-xx[9]][3-xx[8]][3-xx[7]][3-xx[6]][3-xx[5]]*lambda_context_product4[i][tt1[k]][3-xx[4]][3-xx[3]][3-xx[2]][3-xx[1]][3-xx[0]];
							}
							
							//for (int k=0;k<p.length;k++){
							//	System.out.println(p[k]+"	"+lambda[i][k]);//VALID
							//}
						}
						*/
						
						
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					
					//for (int i=0;i<lambda.length;i++){//entities[i]+"	"+
					//	System.out.println(lambda[i][0]+"	"+lambda[i][1]+"	"+lambda[i][2]);
					//}
//					for (int i=0;i<weights.length;i++){
//						for (int j=0;j<weights[i].length;j++){
//							System.out.print("	"+weights[i][j]);
//						}
//						System.out.println();
//					}
					
					//for (int i=0;i<ll.length;i++){//entities[i]+"	"+
					//	System.out.println(ll[i][0]+"	"+ll[i][1]+"	"+ll[i][2]);
					//}
					//System.exit(0);
					
					
//					output.write(pos.get(50)+"	"+nucl.get(50));
//					for (int k=0;k<lambda.length;k++){
//						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
//					}
//					output.newLine();
					
				}
				else if(nucl.get(50).equals("C")){
					//double[][] lambda=new double[lambda_context.size()][3];
					for (int i=0;i<lambda.length;i++){
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][0];
						lambda[i][1]*=lambda_type.get(i)[0][0];
						lambda[i][2]*=lambda_type.get(i)[0][0];
						//int[] tt={0,1,2};
						
						//lambda[i][0]=1;
						//lambda[i][1]=1;
						//lambda[i][2]=1;
						
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
						
						
						/*
						int[] xx=new int[20];
						boolean valid=true;
						for (int j=-1;j>=-offset_left;j--){
							if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
								xx[10+j]=nucl_index(nucl.get(50+j));
							}
							else{
								valid=false;
								break;
							}
						}
						
						if(valid){
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									xx[9+j]=nucl_index(nucl.get(50+j));
								}
								else{
									valid=false;
									break;
								}
							}
						}
						
						if (valid){
							double[] p=new double[3];
							for (int k=0;k<tt2.length;k++){
								
									p[k]=lambda_context_product1[i][tt2[k]][xx[0]][xx[1]][xx[2]][xx[3]][xx[4]]*lambda_context_product2[i][tt2[k]][xx[5]][xx[6]][xx[7]][xx[8]][xx[9]]
										*lambda_context_product3[i][tt2[k]][xx[10]][xx[11]][xx[12]][xx[13]][xx[14]]*lambda_context_product4[i][tt2[k]][xx[15]][xx[16]][xx[17]][xx[18]][xx[19]];
							}
							
							for (int k=0;k<p.length;k++){
	//							System.out.println(p[k]+"	"+lambda[i][k]);//VALID
							}
						}
						*/
						
						
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					
//					lambda2.add(product(weights,lambda));
//					pos2.add(pos.get(50));
//					nucl2.add(nucl.get(50));
//					coverage2.add(coverage.get(50));
//					label2.add(label.get(50));
//					count2.add(count.get(50));
					
//					output.write(pos.get(50)+"	"+nucl.get(50));
//					for (int k=0;k<lambda.length;k++){
//						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
//					}
//					output.newLine();
					
				}
				else if(nucl.get(50).equals("G")){
					//double[][] lambda=new double[lambda_context.size()][3];
					for (int i=0;i<lambda.length;i++){
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][0];
						lambda[i][1]*=lambda_type.get(i)[0][0];
						lambda[i][2]*=lambda_type.get(i)[0][0];
						//int[] tt={0,1,2};
						
						
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
						
						/*
						int[] xx=new int[20];
						boolean valid=true;
						for (int j=-1;j>=-offset_left;j--){
							if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
								xx[10+j]=nucl_index(nucl.get(50+j));
							}
							else{
								valid=false;
								break;
							}
						}
						
						if(valid){
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									xx[9+j]=nucl_index(nucl.get(50+j));
								}
								else{
									valid=false;
									break;
								}
							}
						}
						
						if (valid){
							double[] p=new double[3];
							for (int k=0;k<tt2.length;k++){
								
									p[k]=lambda_context_product1[i][tt2[k]][3-xx[19]][3-xx[18]][3-xx[17]][3-xx[16]][3-xx[15]]*lambda_context_product2[i][tt2[k]][3-xx[14]][3-xx[13]][3-xx[12]][3-xx[11]][3-xx[10]]
										*lambda_context_product3[i][tt2[k]][3-xx[9]][3-xx[8]][3-xx[7]][3-xx[6]][3-xx[5]]*lambda_context_product4[i][tt2[k]][3-xx[4]][3-xx[3]][3-xx[2]][3-xx[1]][3-xx[0]];
							}
							
							for (int k=0;k<p.length;k++){
	//							System.out.println(p[k]+"	"+lambda[i][k]);//VALID
							}
						}
						*/
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					
//					lambda2.add(product(weights,lambda));
//					pos2.add(pos.get(50));
//					nucl2.add(nucl.get(50));
//					coverage2.add(coverage.get(50));
//					label2.add(label.get(50));
//					count2.add(count.get(50));
					
//					output.write(pos.get(50)+"	"+nucl.get(50));
//					for (int k=0;k<lambda.length;k++){
//						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
//					}
//					output.newLine();
					
				}
				else if(nucl.get(50).equals("T")){
					//double[][] lambda=new double[lambda_context.size()][3];
					for (int i=0;i<lambda.length;i++){
						lambda[i][0]=lambda_type.get(i)[1][0];
						lambda[i][1]=lambda_type.get(i)[1][1];
						lambda[i][2]=lambda_type.get(i)[1][2];
						
						lambda[i][0]*=lambda_type.get(i)[0][1];
						lambda[i][1]*=lambda_type.get(i)[0][1];
						lambda[i][2]*=lambda_type.get(i)[0][1];
						//int[] tt={3,5,4};
						
						
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
						
						
						/*
						int[] xx=new int[20];
						boolean valid=true;
						for (int j=-1;j>=-offset_left;j--){
							if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
								xx[10+j]=nucl_index(nucl.get(50+j));
							}
							else{
								valid=false;
								break;
							}
						}
						
						if(valid){
							for (int j=1;j<=offset_right;j++){
								if(pos.get(50)+j==pos.get(50+j)&&nucl_index(nucl.get(50+j))!=-1){
									xx[9+j]=nucl_index(nucl.get(50+j));
								}
								else{
									valid=false;
									break;
								}
							}
						}
						
						if (valid){
							double[] p=new double[3];
							for (int k=0;k<tt2.length;k++){
								
									p[k]=lambda_context_product1[i][tt1[k]][xx[0]][xx[1]][xx[2]][xx[3]][xx[4]]*lambda_context_product2[i][tt1[k]][xx[5]][xx[6]][xx[7]][xx[8]][xx[9]]
										*lambda_context_product3[i][tt1[k]][xx[10]][xx[11]][xx[12]][xx[13]][xx[14]]*lambda_context_product4[i][tt1[k]][xx[15]][xx[16]][xx[17]][xx[18]][xx[19]];
							}
							
							for (int k=0;k<p.length;k++){
//								System.out.println(p[k]+"	"+lambda[i][k]);//VALID
							}
						}
						*/
						
					}
					for (int k=0;k<lambda.length;k++){
						for (int l=0;l<lambda[k].length;l++){
							if(lambda[k][l]>1000){
								lambda[k][l]=1000;
							}
						}
					}
					
					
					
//					lambda2.add(product(weights,lambda));
//					pos2.add(pos.get(50));
//					nucl2.add(nucl.get(50));
//					coverage2.add(coverage.get(50));
//					label2.add(label.get(50));
//					count2.add(count.get(50));
					
//					output.write(pos.get(50)+"	"+nucl.get(50));
//					for (int k=0;k<lambda.length;k++){
//						output.write("	"+lambda[k][0]+";"+lambda[k][1]+";"+lambda[k][2]);
//					}
//					output.newLine();
					
				}
				//System.out.println(x);
				
				
				/*
				lambda2.add(product(weights,lambda));
				pos2.add(pos.get(50));
				nucl2.add(nucl.get(50));
				coverage2.add(coverage.get(50));
				label2.add(label.get(50));
				count2.add(count.get(50));
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
		
		/*
		public static int nucl_index(String s){
			if(s.toUpperCase().equals("A")){
				return 0;
			}
			else if(s.toUpperCase().equals("C")){
				return 1;
			}
			else if(s.toUpperCase().equals("G")){
				return 2;
			}
			else if(s.toUpperCase().equals("T")){
				return 3;
			}
			return -1;
		}*/
		
		public static int nucl_index(String s){
			Integer ii=table_nucl.get(s);
			if(ii!=null){
				return ii.intValue();
			}
			else{
				return -1;
			}
		}
	
	
	public static double[][] product(double[][] x, double[][] y){
		double[][] z=new double[x.length][y[0].length];
		for (int i=0;i<x.length;i++){
			for (int j=0;j<y[0].length;j++){
				for (int k=0;k<y.length;k++){
					z[i][j]+=x[i][k]*y[k][j];
				}
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
}
