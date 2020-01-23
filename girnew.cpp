#include<bits/stdc++.h>
using namespace std;

class Graph{ 
    //int V; //nodes
    list<pair<int,int> >* adj; 
    vector<int>deg;
    vector<pair<int,int> > edges;
    void dfs(int v, bool visited[]); 
    bool weight;
  	
public: 	
	int V;
    Graph(int V); 
  	double m_;
  	int edgesno;
    void addEdge(int u,int v, int w); 
    int NumberOfconnectedComponents(); 
    void gnstep();
    void remove_edge(int u,int v);
    map<pair<int,int>,double> edgebc();
    void cc(int v,bool visited[],double &ewc,double &re,vector<int>deg1);
    double gnMod();
    vector<int> update();
    void rungn();
    void bfspath(vector<int> &S,vector<int> P[],map<int,double> &sigma,int s);
    void dijkspath(vector<int> &S,vector<int> P[],map<int,double> &sigma,int s);
    void acc_edges(map<pair<int,int>,double> &bet,vector<int>&S,vector<int>P[],map<int,double>&sigma,int s);
    void rescale_e(map<pair<int,int>,double>&bet);
    void answerutil(int v,vector<bool>&vis,vector<int>&comp);
    void answer(vector<vector<int> > &ans);
}; 
   

void Graph::remove_edge(int u,int v){
	int e=edges.size();
	for(int i=e-1;i>=0;--i){
		if((edges[i].first==u && edges[i].second==v) || (edges[i].first==v && edges[i].second==u)){
			edges.erase(edges.begin()+i);
			edgesno--;
			continue;
		}
	}
	list<pair<int,int> >a=adj[v];
	for(auto it=a.begin();it!=a.end();it++){
		if(it->first==u){
			it = a.erase(it);
			continue;
		}
	}
	adj[v]=a;
	a=adj[u];
	for(auto it=a.begin();it!=a.end();it++){
		if(it->first==v){
			it=a.erase(it);
			continue;
		}
	}
	adj[u]=a;
}

void Graph::gnstep(){
    int init_comp=NumberOfconnectedComponents();
    int ncomp=init_comp;
    while(ncomp<=init_comp){
    	//	cout<<ncomp<<" "<<init_comp<<endl;
        map<pair<int,int>,double> bet=edgebc();
        //cout<<bet.size()<<endl;
        double ma=-1000000.0;
        for(auto it=bet.begin();it!=bet.end();it++){
        	pair<int,int> pa=it->first;
        	//cout<<pa.first<<" "<<pa.second<<" "<<it->second<<endl;
            if(ma<it->second){
                ma=it->second;
            }
        }
        //cout<<ma<<endl;
        for(auto it=bet.begin();it!=bet.end();it++){
            if(it->second==ma){
                pair<int,int> pt=it->first;
                //cout<<pt.first<<" "<<pt.second<<" "<<ma<<endl;
                remove_edge(pt.first,pt.second);
                //cout<<"lopop"<<endl;
            }
        }
        //cout<<"ll"<<endl;
        ncomp = NumberOfconnectedComponents();
        //cout<<ncomp<<" "<<init_comp<<endl;
    }
}


void Graph::cc(int v,bool visited[],double &ewc,double &re,vector<int>deg1){
	visited[v]=true;
	ewc+=deg1[v];
	re+=deg[v];
	for(auto it=adj[v].begin();it!=adj[v].end();it++){
		pair<int,int> p=*it;
		if(!visited[p.first]){
			cc(p.first,visited,ewc,re,deg1);
		}
	}
}

double Graph::gnMod(){
	vector<int>deg1 = update();
	int comps = NumberOfconnectedComponents();
	cout<<"Number of connectedcomponents: "<<comps<<endl;
	double mod=0.0;
	bool *visited = new bool[V];
	for(int i=0;i<V;i++){
		visited[i]=false;
	}
	for(int i=0;i<V;i++){
		double ewc=0.0;
		double re=0.0;
		if(visited[i]==false){
			cc(i,visited,ewc,re,deg1);
		}
		mod+=(ewc)-((re*re)*1.0)/(2.0*m_);
	}
	mod=mod/(2.0*m_);
	return mod;
}

void Graph::answerutil(int v,vector<bool>&vis,vector<int>&comp){
	vis[v]=true;
	comp.push_back(v);
	for(auto i=adj[v].begin();i!=adj[v].end();i++){
		if(!vis[i->first]){
			answerutil(i->first,vis,comp);
		}
	}
}


void Graph::answer(vector<vector<int> > &ans){
	vector<bool> vis(V,false);
	vector<int> comp;
	for(int i=0;i<V;i++){
		if(vis[i]==false){
			answerutil(i,vis,comp);
			if(comp.size()!=0){
				ans.push_back(comp);
			}
			comp.clear();
		}
	}
}


void Graph::rungn(){
	double bestQ=0.0;
	double Q=0.0;
	vector<vector<int> > ans;
	while(1){
		gnstep();
		Q=gnMod();
		cout<<"Modularity of Decomposed G: "<<Q<<endl;
		if(Q>bestQ){
			ans.clear();
			bestQ=Q;
			answer(ans);
			for(int i=0;i<ans.size();i++){
				for(int j=0;j<ans[i].size();j++){
					cout<<ans[i][j]<<" ";
				}
				cout<<endl;
			}
			//cout<<"answer"<<endl;
		}
		if(edgesno==0){
			break;
		}
		//cout<<edgesno<<endl;
	}
	if(bestQ>0.0){
		//answeru
		cout<<"Max modularity(Q): "<<bestQ<<endl;
		for(int i=0;i<ans.size();i++){
			for(int j=0;j<ans[i].size();j++){
				cout<<ans[i][j]<<" ";
			}
			cout<<endl;
		}
		//cout<<"answer"<<endl;
	}
	else{
		cout<<"Modularity of Decomposed G: "<<Q<<endl;
	}
}

vector<int> Graph::update(){
	vector<int>temp(V,0);
	for(int i=0;i<V;i++){
		list<pair<int,int> > a=adj[i];
		for(auto j=a.begin();j!=a.end();j++){
			temp[i]+=j->second;
		}
	}
	return temp;
}

map<pair<int,int>,double> Graph::edgebc(){
    map<pair<int,int>,double> bet;
    double ze=0.0;
    /*for(int i=0;i<V;i++){
        bet[{i,-1}]=ze;
    }*/
    for(int j=0;j<edges.size();j++){
        bet[{edges[j].first,edges[j].second}]=ze;
    }
    for(int s=0;s<V;s++){
        vector<int> S;
        vector<int>P[V];
        map<int,double>sigma;
        for(int i=0;i<V;i++){
        	sigma[i]=0.0;
        }
        if(weight){
            bfspath(S,P,sigma,s);
        }
        else{
            dijkspath(S,P,sigma,s);
        }
        acc_edges(bet,S,P,sigma,s);
    }
    /*for(int i=0;i<V;i++){
    	auto it = bet.find({i,-1});
    	if(it!=bet.end()){
    		bet.erase(it);	
    	}
    }*/
    //rescale_e(bet);
    return bet;
}


void Graph::acc_edges(map<pair<int,int>,double> &bet,vector<int>&S,vector<int>P[],map<int,double>&sigma,int s){
    map<int,double> delta;
    for(int i=0;i<S.size();i++){
        delta[S[i]]=0.0;
    }
    int l=S.size()-1;
    while(l>=0){
        int w=S[l--];
        double coeff=(1.0+delta[w])/(sigma[w]*1.0);
        for(int i=0;i<P[w].size();i++){
            int v=P[w][i];
            double c=sigma[v]*coeff;
            if(bet.find({v,w})!=bet.end()){
                bet[{v,w}]+=c;
            }
            else{
                bet[{w,v}]+=c;
            }
            delta[v]+=c;
        }
        /*if(w!=s){
            bet[{w,-1}]+=delta[w];
        }*/
    }
}

void Graph::rescale_e(map<pair<int,int>,double>&bet){
    double scale;
    double l=1.0;
    scale=l/(V*(V-l));
    
    for(auto it=bet.begin();it!=bet.end();it++){
        pair<int,int> pt=it->first;
    	bet[pt]*=scale;
	}
}

void Graph::bfspath(vector<int> &S,vector<int> P[],map<int,double> &sigma,int s){
    map<int,int> D;
    double on=1.0;
    sigma[s]=1.0;
    D[s]=0;
    queue<int> q;
    q.push(s);
    while(!q.empty()){
        int v=q.front();
        q.pop();
        S.push_back(v);
        int Dv=D[v];
        double sigmav=sigma[v];
        list<pair<int,int> > a=adj[v];
        for(auto i=a.begin();i!=a.end();i++){
            int w=i->first;
            if(D.find(w)==D.end()){
                q.push(w);
                D[w]=Dv+1;
            }
            if(D[w]==Dv+1){
                sigma[w]+=sigmav;
                P[w].push_back(v);
            }
        }
    }
}


void Graph::dijkspath(vector<int> &S,vector<int> P[],map<int,double> &sigma,int s){
    map<int,int> seen;
    map<int,int> D;
    seen[s]=0;
    int c=0;
    sigma[s]=1.0;
    queue<vector<int> > q;
    vector<int> v;
    v.push_back(0);
    v.push_back(c++);
    v.push_back(s);
    v.push_back(s);
    q.push(v);
    while(!q.empty()){
        vector<int> temp=q.front();
        q.pop();
        int dist=temp[0];
        int pred=temp[2];
        int v=temp[3];
        if(D.find(v)!=D.end()){
            continue;
        }
        sigma[v]+=sigma[pred];
        S.push_back(v);
        D[v]=dist;
        list<pair<int,int> > adjl=adj[v];
        for(auto i=adjl.begin();i!=adjl.end();i++){
            int w=i->first;
            int wt=i->second;
            int vw_dist=dist+wt;
            if(D.find(w)==D.end() && (seen.find(w)==seen.end() || vw_dist < seen[w])){
                seen[w]=vw_dist;
                vector<int> vt;
                vt.push_back(vw_dist);
                vt.push_back(c++);
                vt.push_back(v);
                vt.push_back(w);
                q.push(vt);
                sigma[w]=0.0;
                P[w].clear();
                P[w].push_back(v);
            }
            else if(vw_dist==seen[w]){
                sigma[w]+=sigma[v];
                P[w].push_back(v);
            }
        }
    }
}


int Graph::NumberOfconnectedComponents(){ 
    bool* visited = new bool[V];
    int count = 0; 
    for (int v = 0; v < V; v++) 
        visited[v] = false; 
  
    for (int v = 0; v < V; v++){ 
        if (visited[v] == false){
            dfs(v, visited); 
            count += 1; 
        } 
    } 
  
    return count; 
} 
  
void Graph::dfs(int v, bool visited[]){ 
    visited[v] = true;
    for (auto i = adj[v].begin(); i != adj[v].end(); ++i) 
        if (!visited[i->first]) 
            dfs(i->first, visited); 
}
  
Graph::Graph(int V){ 
    this->V = V; 
    adj = new list<pair<int,int> >[V]; 
    deg.resize(V);
    m_=0;
    edgesno=0;
    weight=1;
} 

void Graph::addEdge(int u,int v, int w){ 
    adj[v].push_back({u,w}); 
    adj[u].push_back({v,w}); 
    edges.push_back({u,v});
    deg[u]+=w;
    deg[v]+=w;
    m_+=w;
    edgesno++;
} 

int main(){
    int n,e;
    cin>>n>>e;
    Graph g(n);
    for(int i=0;i<e;i++){
    	int x,y;
    	cin>>x>>y;
    	g.addEdge(x,y,1);
    }
    g.rungn();
    return 0;
}
