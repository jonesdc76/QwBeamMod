{
  for(int i=44; i<137;++i){
    TFile *f = TFile::Open(Form("excluded_quartets/reduced_slug%i_pattnum.root",i));
    if(f == 0) continue;
    cout<<i<<" "<<slug->GetEntries()<<endl;
  }


}
