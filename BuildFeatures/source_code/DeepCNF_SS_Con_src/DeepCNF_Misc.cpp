#include "DeepCNF_Misc.h"

// ---- parse string -----//
int Parse_Str(string &in,vector <int> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		int value=atoi(buf.c_str());
		out.push_back(value);
		count++;
	}
	return count;
}
int Parse_Str_Double(string &in,vector <double> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		double value=atof(buf.c_str());
		out.push_back(value);
		count++;
	}
	return count;
}
void Parse_Double(vector <double> &in, string &out,char separator)
{
	stringstream oss;
	string sepa="";
	sepa=sepa+separator;
	for(int i=0;i<(int)in.size();i++)
	{
		int wsiii=(int)in[i];
		if(wsiii!=in[i])oss << in[i] << sepa;
		else oss << wsiii << sepa;
	}
	string wsbuf=oss.str();
	out=wsbuf;
}

//------- parse feature range --------//
//-> string to int
int str2int_(const char *str,int &num)
{
	int i,j,len=(int)strlen(str);
	for(i=0;i<len;i++)if(str[i]!=' ')break;
	for(j=len-1;j>=i;j--)if(str[j]!=' ')break;
	if(str[i]=='-')i++;
	for(;i<=j;i++)if(str[i]<'0'||str[i]>'9')return 0; //error
	num=atoi(str);
	return 1; //correct
}

//-> parse a single range, such as 1-7, or 9
void Parse_Single_Range(string &in_str,int &start,int &end)
{
	int pos=-1;
	for(int i=0;i<(int)in_str.length();i++)
	{
		if(in_str[i]=='-')
		{
			pos=i;
			break;
		}
	}
	if(pos==-1) //-> single position
	{
		int value;
		int retv=str2int_(in_str.c_str(),value);
		if(retv!=1)
		{
			fprintf(stderr,"range %s error !! \n",in_str.c_str());
			exit(-1);
		}
		start=value-1;
		end=value-1;
	}
	else        //-> start to end
	{
		string start_str=in_str.substr(0,pos);
		string end_str=in_str.substr(pos+1,in_str.length()-pos-1);
		int retv;
		int value;
		//--> get start
		retv=str2int_(start_str.c_str(),value);
		if(retv!=1)
		{
			fprintf(stderr,"range start part %s error !! \n",start_str.c_str());
			exit(-1);
		}
		start=value-1;
		//--> get end
		retv=str2int_(end_str.c_str(),value);
		if(retv!=1)
		{
			fprintf(stderr,"range end part %s error !! \n",end_str.c_str());
			exit(-1);
		}
		end=value-1;
	}
}

//-> we allow the feature selection by the following format: 1-7,9,24-30 //starting from 1
int Parse_Feature_Range(string &in_str_,vector <int> &range_out)
{
	string in_str=in_str_+",";
	for(int i=0;i<(int)range_out.size();i++)range_out[i]=0;
	//get block
	vector <string> init_block;
	istringstream www(in_str);
	for(;;)
	{
		string buf;
		if(!getline(www,buf,','))break;
		init_block.push_back(buf);
	}
	//analyse block
	for(int i=0;i<(int)init_block.size();i++)
	{
		//-> get start and end
		int start,end;
		Parse_Single_Range(init_block[i],start,end);
		if(start<0 || end>=range_out.size())
		{
			fprintf(stderr,"WARNING: range_string %s overrange : start %d, end %d, length %d !! \n",
				init_block[i].c_str(),start,end,(int)range_out.size());
			//--> fix the problem
			if(start<0)start=0;
			if(end>=range_out.size())end=range_out.size()-1;
		}
		//-> assign value
		for(int k=start;k<=end;k++)range_out[k]=1;
	}
	//count remaining feature number
	int feat_num=0;
	for(int i=0;i<(int)range_out.size();i++)if(range_out[i]==1)feat_num++;
	return feat_num;
}

//==================== load training data =====================//
//-> if with single CPU, set num_procs=1 and proc_id=0
//-> data format ( length, features, labels ) 
/*
78
feat1 feat2 feat3, ..., featN
....  (with 78 feature lines)
0
1
0
2
...  ( with 78 label lines)
54
feat1 feat2 feat3, ..., featN
....  (with 54 feature lines)
0
1
0
2
...  ( with 54 label lines)
*/

//---- parse feat_num ---//
int Parse_FeatNum(string &in_str)
{
	int feat_num=0;
	istringstream www(in_str);
	for(;;)
	{
		string temp;
		if(! (www>>temp) )break;
		feat_num++;
	}
	return feat_num;
}

//---- parse selection ---//
void Parse_Selection(string &in_str, string &out_str, vector <int> &range_out)
{
	//inic check
	out_str=in_str;
	int size=(int)range_out.size();
	int rel=0;
	for(int i=0;i<size;i++)if(range_out[i]==1)rel++;
	if(rel==size)return;

	//load feature
	istringstream www(in_str);
	vector <string> tmp_rec;
	for(;;)
	{
		string temp;
		if(! (www>>temp) )break;
		tmp_rec.push_back(temp);
	}
	//get new feature
	out_str="";
	for(int i=0;i<(int)range_out.size();i++)
	{
		if(range_out[i]==1)
		{
			out_str=out_str+tmp_rec[i]+" ";
		}
	}
}

//---- load data ----//
int LoadData(string &input_file, int num_procs,int proc_id, 
	int local_num_ori, vector <int> &range_out,
	vector <vector <string> > & feat_in, vector <vector <int> > & label_in)
{
	ifstream fin;
	string buf,temp;
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	//read in total data
	feat_in.clear();
	label_in.clear();
	int count=-1;
	int feat_num_ori=local_num_ori;
	for(;;)
	{
		//-> construct a new sequence
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		int length_s=atoi(buf.c_str());
		count++;
		//-> read in features
		vector <string> feat_rec;
		for(int j=0;j<length_s;j++)
		{
			if(!getline(fin,buf,'\n'))
			{
				fprintf(stderr,"input_file %s file format bad!\n",input_file.c_str());
				exit(-1);
			}
			//--> check feature number
			int feat_num=Parse_FeatNum(buf);
			if(feat_num!=feat_num_ori)
			{
				fprintf(stderr,"current feat_num %d not equal to first feat_num %d at line %s \n",
					feat_num,feat_num_ori,buf.c_str());
				exit(-1);
			}
			string buf_;
			Parse_Selection(buf, buf_, range_out);
			//--> push_back
			feat_rec.push_back(buf_);
		}
		//-> read in labels
		vector <int> label_rec;
		for(int j=0;j<length_s;j++)
		{
			if(!getline(fin,buf,'\n'))
			{
				fprintf(stderr,"input_file %s file format bad!\n",input_file.c_str());
				exit(-1);
			}
			label_rec.push_back(atoi(buf.c_str()));
		}
		//-> for MPI purpose
		if(count%num_procs != proc_id)continue;
		//-> add to data
		feat_in.push_back(feat_rec);
		label_in.push_back(label_rec);
	}
	//return count
	return count+1;
}
int LoadData(string &input_file, int num_procs,int proc_id, 
	int local_num_ori, vector <int> &range_out,
	vector <vector <vector <FeatScore> > > & feat_in, vector <vector <int> > & label_in)
{
	ifstream fin;
	string buf,temp;
	fin.open(input_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"input_file %s not found!\n",input_file.c_str());
		exit(-1);
	}
	//read in total data
	feat_in.clear();
	label_in.clear();
	int count=-1;
	int feat_num_ori=local_num_ori;
	for(;;)
	{
		//-> construct a new sequence
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		int length_s=atoi(buf.c_str());
		count++;
		//-> read in features
		vector <vector <FeatScore> > feat_rec;
		for(int j=0;j<length_s;j++)
		{
			if(!getline(fin,buf,'\n'))
			{
				fprintf(stderr,"input_file %s file format bad!\n",input_file.c_str());
				exit(-1);
			}
			//--> check feature number
			int feat_num=Parse_FeatNum(buf);
			if(feat_num!=feat_num_ori)
			{
				fprintf(stderr,"current feat_num %d not equal to first feat_num %d at line %s \n",
					feat_num,feat_num_ori,buf.c_str());
				exit(-1);
			}
			string buf_;
			Parse_Selection(buf, buf_, range_out);
			//--> process data
			vector <FeatScore> value_rec;
			istringstream trn_in(buf_);
			for(;;)
			{
				FeatScore value;
				if(! (trn_in >> value) )break;
				value_rec.push_back(value);
			}
			feat_rec.push_back(value_rec);
		}
		//-> read in labels
		vector <int> label_rec;
		for(int j=0;j<length_s;j++)
		{
			if(!getline(fin,buf,'\n'))
			{
				fprintf(stderr,"input_file %s file format bad!\n",input_file.c_str());
				exit(-1);
			}
			label_rec.push_back(atoi(buf.c_str()));
		}
		//-> for MPI purpose
		if(count%num_procs != proc_id)continue;
		//-> add to data
		feat_in.push_back(feat_rec);
		label_in.push_back(label_rec);
	}
	//return count
	return count+1;
}
