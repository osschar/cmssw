#include "Fireworks/Core/src/FWEveWriteToObjFile.h"
#include "Fireworks/Core/interface/FWEveView.h"

#include "TEveCalo.h"
#include "TEveCaloData.h"
#include "TEveGeoShape.h"
#include "TEveProjections.h"
#include "TEveProjectionManager.h"
#include "TEveRGBAPalette.h"
#include "TEveText.h"
#include "TEveTrans.h"
#include "TEveTrack.h"

#include "TGLUtil.h"
#include "TGListTree.h"

#include "TGeoBBox.h"
#include "TClass.h"
#include "TMathBase.h"
#include "TMath.h"
#include "TAxis.h"

#include <iostream>
#include <fstream>
#include <stdexcept>

FWEveWriteToObjFile::FWEveWriteToObjFile(TGListTreeItem * first_item_of_full_tree)
{
  fM=0;
  tracks_from_3dtower=0;
  ecal_from_3dtower=0;
  muons_from_3dtower=0;
  transF_eta = 1.479;
  transB_eta = -1.479;
  endcapF_z = 0;
  endcapB_z = -0; 

  scale=1.0;

  n_vert=0;
  n_norm=0;

  eventscene_3dtower = my_FindSiblingByName(first_item_of_full_tree,"EventScene 3D Tower");
  if (!eventscene_3dtower){
    std::cout<<"Can't find \"EventScene 3D Tower\" !!!!"<<std::endl;
    return;
  }
  
  tracks_from_3dtower = my_FindSiblingByName(eventscene_3dtower,"TracksFromPFCands");
  if (!tracks_from_3dtower)
    tracks_from_3dtower=my_FindSiblingByName(eventscene_3dtower,"Tracks");

  ecal_from_3dtower=my_FindSiblingByName(eventscene_3dtower,"calo barrel");
  if (ecal_from_3dtower) {
    fM =(TEveCalo3D*) ecal_from_3dtower->GetUserData();
    //For CMSSW_5_X_Y
    //transF_eta = fM->GetTransitionEta();
    //transB_eta = -fM->GetTransitionEta();
    //endcapF_z = fM->fEndCapPos;
    //endcapB_z = -fM->fEndCapPos; 

    // For CMSSW_7_X_Y
    transF_eta  = fM->GetTransitionEtaForward();
    transB_eta  = fM->GetTransitionEtaBackward();
    endcapF_z = fM->GetForwardEndCapPos();
    endcapB_z = fM->GetBackwardEndCapPos();
  }
  muons_from_3dtower=my_FindSiblingByName(eventscene_3dtower,"Muons");
}


void FWEveWriteToObjFile::save_to_obj(const char* iName ,bool save_tracks, bool save_calo, bool save_muon_track, bool save_muon_chamber)
{
  n_vert=0;
  n_norm=0;
  std::ofstream out_file(iName);
  if (save_tracks){
    if (!tracks_from_3dtower) std::cout<<"tracks are not found"<<std::endl;
    else{
      int written_tracks=1;
      write_tracks(out_file,"track_",tracks_from_3dtower->GetFirstChild(),written_tracks);
    }
  }
  if (save_calo){
    if (!fM) std::cout<<"Calo information are not found"<<std::endl;
    else write_calo(out_file,"Calo_");
  }

  if (save_muon_track){
    if (!muons_from_3dtower) std::cout<<"muon tracks are not found"<<std::endl;
    else {
      int written_tracks=1;
      write_tracks(out_file,"muon_track_",muons_from_3dtower->GetFirstChild(),written_tracks);    
    }
  }
  if (save_muon_chamber){
    if (!muons_from_3dtower) std::cout<<"muon chambers are not found"<<std::endl;
    else {
      int written_chambers=1;
      write_chambers(out_file,"muon_chamb_",muons_from_3dtower->GetFirstChild(),written_chambers);    
    }
  }
  out_file.close();
}

TGListTreeItem * FWEveWriteToObjFile::my_FindSiblingByName(TGListTreeItem *item, const TString& name)
{
  while (item) {
    TString akt_name=item->GetText();
    if (akt_name.EqualTo(name)) return item;
    if (item->GetFirstChild()!=0) {
      TGListTreeItem * test_item=my_FindSiblingByName(item->GetFirstChild(),name);
      if (test_item !=0) return test_item;
    }
    item=item->GetNextSibling();
  }
  return 0;
}

//==============================================================================

TString FWEveWriteToObjFile::give_track_vertecies(TEveTrack* akt_track)
{
  if (!akt_track) return " ";
  std::stringstream  ret_str;
  Float_t xx,yy,zz;
  for (int i=0; i<akt_track->GetN(); i++){
    akt_track->GetPoint(i,xx,yy,zz);
    ret_str<<"v "<<xx*scale<<"\t"<<yy*scale<<"\t"<<zz*scale<<"\n";
    n_vert++;
  }
  return ret_str.str();
}


void FWEveWriteToObjFile::write_tracks(std::ofstream &out_file, const TString& track_name, TGListTreeItem *akt_tree_item, int &n_track)
{
  while (akt_tree_item){
    if (TString(akt_tree_item->GetText()).Contains("TEveRecTrackT")){
      TEveTrack* akt_track=(TEveTrack*) akt_tree_item->GetUserData();
      out_file<<"o "<<track_name<<n_track<<'\n';
      out_file<<give_track_vertecies(akt_track);
      n_track++;
    }
    if (akt_tree_item->GetFirstChild()!=0) write_tracks(out_file,track_name,akt_tree_item->GetFirstChild(),n_track);
    akt_tree_item=akt_tree_item->GetNextSibling();
  }
}


void FWEveWriteToObjFile::write_chambers(std::ofstream &out_file, const TString& chamber_name, TGListTreeItem *akt_tree_item, int &n_chamber)
{
  while (akt_tree_item){
    if (TString(akt_tree_item->GetText()).Contains("Chamber")){
      TEveGeoShape* akt_shape=(TEveGeoShape*) akt_tree_item->GetUserData();
      TGeoBBox*  akt_geo_box=(TGeoBBox*) akt_shape->GetShape();

      out_file<<"o "<<chamber_name<<n_chamber<<'\n';
      TEveTrans& trans_form=akt_shape->RefMainTrans();
      double local_box[24];
      akt_geo_box->SetPoints(local_box);
      double* local_pbox=local_box;
      float global_box[24];
      for (int ii=0; ii<8; ii++){
	trans_form.MultiplyIP(local_pbox);
	global_box[3*ii]=*local_pbox;
	global_box[3*ii+1]=*(local_pbox+1);
	global_box[3*ii+2]=*(local_pbox+2);
	local_pbox+=3;
      }
      out_file<<give_box_vetrecies_and_norm(global_box);
      n_chamber++;
    }
    if (akt_tree_item->GetFirstChild()!=0) write_chambers(out_file,chamber_name,akt_tree_item->GetFirstChild(),n_chamber);
    akt_tree_item=akt_tree_item->GetNextSibling();
  }
}



TString FWEveWriteToObjFile::give_box_vetrecies_and_norm(Float_t box[24]){
 
  std::stringstream  ret_str;
  //    z
  //    |
  //    |
  //    |________y
  //   /  6-------7
  //  /  /|      /|
  // x  5-------4 |
  //    | 2-----|-3
  //    |/      |/
  //    1-------0
  //

  for (int ii=0; ii<8; ii++){
    ret_str<<"v "<<box[3*ii]*scale<<"\t"<<box[3*ii+1]*scale<<"\t"<<box[3*ii+2]*scale<<"\n";
    n_vert++;
  }


  Float_t cross[3];
  const Float_t *p = box;
  double abs_cross=1;
  double norm_scale=1;


  // bottom: 0123
  CrossProduct(p+3, p+9, p, cross);
  abs_cross=(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  if (abs_cross<=0) norm_scale=0;
  else norm_scale=1/sqrt(abs_cross);
  ret_str<<"vn "<<double(cross[0])*norm_scale<<"\t"<<double(cross[1])*norm_scale<<"\t"<<double(cross[2])*norm_scale<<"\n";
  n_norm++;
  // top:    7654
  CrossProduct(p+21, p+15, p+12, cross);
  abs_cross=(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  if (abs_cross<=0) norm_scale=0;
  else norm_scale=1/sqrt(abs_cross);
  ret_str<<"vn "<<double(cross[0])*norm_scale<<"\t"<<double(cross[1])*norm_scale<<"\t"<<double(cross[2])*norm_scale<<"\n";
  n_norm++;
  // back:   0451
  CrossProduct(p+12, p+3, p, cross);
  abs_cross=(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  if (abs_cross<=0) norm_scale=0;
  else norm_scale=1/sqrt(abs_cross);
  ret_str<<"vn "<<double(cross[0])*norm_scale<<"\t"<<double(cross[1])*norm_scale<<"\t"<<double(cross[2])*norm_scale<<"\n";
  n_norm++;
  //front :  3267
  CrossProduct(p+6, p+21, p+9, cross);
  abs_cross=(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  if (abs_cross<=0) norm_scale=0;
  else norm_scale=1/sqrt(abs_cross);
  ret_str<<"vn "<<double(cross[0])*norm_scale<<"\t"<<double(cross[1])*norm_scale<<"\t"<<double(cross[2])*norm_scale<<"\n";
  n_norm++;
  // left:    0374
  CrossProduct(p+21, p, p+9, cross);
  abs_cross=(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  if (abs_cross<=0) norm_scale=0;
  else norm_scale=1/sqrt(abs_cross);
  ret_str<<"vn "<<double(cross[0])*norm_scale<<"\t"<<double(cross[1])*norm_scale<<"\t"<<double(cross[2])*norm_scale<<"\n";
  n_norm++;
  // right:   1562
  CrossProduct(p+15, p+6, p+3, cross);
  abs_cross=(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  if (abs_cross<=0) norm_scale=0;
  else norm_scale=1/sqrt(abs_cross);
  ret_str<<"vn "<<double(cross[0])*norm_scale<<"\t"<<double(cross[1])*norm_scale<<"\t"<<double(cross[2])*norm_scale<<"\n";
  n_norm++;

  int p0=(n_vert-7);
  int p1=(n_vert-6);
  int p2=(n_vert-5);
  int p3=(n_vert-4);
  int p4=(n_vert-3);
  int p5=(n_vert-2);
  int p6=(n_vert-1);
  int p7=(n_vert);

  int nbottom=(n_norm-5);
  int ntop=(n_norm-4);
  int nback=(n_norm-3);
  int nfront=(n_norm-2);
  int nleft=(n_norm-1);
  int nright=(n_norm);

  // face bottom  0123
  ret_str<<"f "<<p0<<"//"<<nbottom<<"\t"<<p1<<"//"<<nbottom<<"\t"<<p2<<"//"<<nbottom<<"\n";
  ret_str<<"f "<<p0<<"//"<<nbottom<<"\t"<<p2<<"//"<<nbottom<<"\t"<<p3<<"//"<<nbottom<<"\n";
  // face top 7654 
  ret_str<<"f "<<p7<<"//"<<ntop<<"\t"<<p6<<"//"<<ntop<<"\t"<<p5<<"//"<<ntop<<"\n";
  ret_str<<"f "<<p7<<"//"<<ntop<<"\t"<<p5<<"//"<<ntop<<"\t"<<p4<<"//"<<ntop<<"\n";
  // face back 0451
  ret_str<<"f "<<p0<<"//"<<nback<<"\t"<<p4<<"//"<<nback<<"\t"<<p5<<"//"<<nback<<"\n";
  ret_str<<"f "<<p0<<"//"<<nback<<"\t"<<p5<<"//"<<nback<<"\t"<<p1<<"//"<<nback<<"\n";
  // face front 3267
  ret_str<<"f "<<p3<<"//"<<nfront<<"\t"<<p2<<"//"<<nfront<<"\t"<<p6<<"//"<<nfront<<"\n";
  ret_str<<"f "<<p3<<"//"<<nfront<<"\t"<<p6<<"//"<<nfront<<"\t"<<p7<<"//"<<nfront<<"\n";
  // face left 0374
  ret_str<<"f "<<p0<<"//"<<nleft<<"\t"<<p3<<"//"<<nleft<<"\t"<<p7<<"//"<<nleft<<"\n";
  ret_str<<"f "<<p0<<"//"<<nleft<<"\t"<<p7<<"//"<<nleft<<"\t"<<p4<<"//"<<nleft<<"\n";
  // face right 1562
  ret_str<<"f "<<p1<<"//"<<nright<<"\t"<<p5<<"//"<<nright<<"\t"<<p6<<"//"<<nright<<"\n";
  ret_str<<"f "<<p1<<"//"<<nright<<"\t"<<p6<<"//"<<nright<<"\t"<<p2<<"//"<<nright<<"\n";

  return ret_str.str();
}





void FWEveWriteToObjFile::write_calo(std::ofstream &out_file, const TString& calo_name) ///  == DirectDraw() from  graf3d/eve/src/TEveCalo3DGL.cxx
{

  Float_t towerH = 0;
  int N_hcal=0;
  int N_ecal=0;
  Int_t   tower = 0;
  Int_t   prevTower = -1;
  Float_t offset = 0;
  Int_t cellID = 0;


  TEveCaloData::vCellId_t akt_cells;
  fM->GetData()->GetCellList(0,fM->GetEtaMax(),fM->GetPhi(),fM->GetPhiRng(),akt_cells);
  fOffset.assign(akt_cells.size(), 0);
  for (TEveCaloData::vCellId_i it=akt_cells.begin(); it!=akt_cells.end(); it++){
    TEveCaloData::CellData_t cellData;
    fM->GetData()->GetCellData((*it),cellData);
    out_file<<"o "<<calo_name;
    tower = it->fTower;
    if (tower != prevTower)
      {
	offset = 0;
	prevTower = tower;
      }
    fOffset[cellID] = offset;

    if (fM->GetValueIsColor())
    {
      towerH = fM->GetValToHeight() * fM->GetData()->GetMaxVal(fM->GetPlotEt());
    }
    else
    {
      towerH = fM->GetValToHeight() * cellData.Value(fM->GetPlotEt());
    }

    if (it->fSlice == 1){
       N_ecal++;
       out_file<<"ECAL_"<<N_ecal; 
    }
    if (it->fSlice == 2){
      N_hcal++; 
      out_file<<"HCAL_"<<N_hcal;
    }
    out_file<<"\n"; 

    Float_t box[24];
    if ((cellData.Eta() > 0 && cellData.Eta() < transF_eta) ||
        (cellData.Eta() < 0 && cellData.Eta() > transB_eta))
    {
      RenderBarrelCell(cellData, towerH, offset,box);
      //	out_file<<give_box_vetrecies_and_norm(box);
    }
    else
    {
      RenderEndCapCell(cellData, towerH, offset,box);
    }  
    out_file<<give_box_vetrecies_and_norm(box);
  }
}

//______________________________________________________________________________
void FWEveWriteToObjFile::RenderBarrelCell(const TEveCaloData::CellGeom_t &cellData, Float_t towerH, Float_t& offset, Float_t* pnts) const // 
{
   // Render barrel cell.

   using namespace TMath;

   Float_t r1 = fM->GetBarrelRadius() + offset;
   Float_t r2 = r1 + towerH*Sin(cellData.ThetaMin());
   Float_t z1In, z1Out, z2In, z2Out;

   z1In  = r1/Tan(cellData.ThetaMax());
   z1Out = r2/Tan(cellData.ThetaMax());
   z2In  = r1/Tan(cellData.ThetaMin());
   z2Out = r2/Tan(cellData.ThetaMin());

   Float_t cos1 = Cos(cellData.PhiMin());
   Float_t sin1 = Sin(cellData.PhiMin());
   Float_t cos2 = Cos(cellData.PhiMax());
   Float_t sin2 = Sin(cellData.PhiMax());

   //   Float_t box[24];
   //Float_t* pnts = box;
   // 0
   pnts[0] = r1*cos2;
   pnts[1] = r1*sin2;
   pnts[2] = z1In;
   pnts += 3;
   // 1
   pnts[0] = r1*cos1;
   pnts[1] = r1*sin1;
   pnts[2] = z1In;
   pnts += 3;
   // 2
   pnts[0] = r1*cos1;
   pnts[1] = r1*sin1;
   pnts[2] = z2In;
   pnts += 3;
   // 3
   pnts[0] = r1*cos2;
   pnts[1] = r1*sin2;
   pnts[2] = z2In;
   pnts += 3;
   //---------------------------------------------------
   // 4
   pnts[0] = r2*cos2;
   pnts[1] = r2*sin2;
   pnts[2] = z1Out;
   pnts += 3;
   // 5
   pnts[0] = r2*cos1;
   pnts[1] = r2*sin1;
   pnts[2] = z1Out;
   pnts += 3;
   // 6
   pnts[0] = r2*cos1;
   pnts[1] = r2*sin1;
   pnts[2] = z2Out;
   pnts += 3;
   // 7
   pnts[0] = r2*cos2;
   pnts[1] = r2*sin2;
   pnts[2] = z2Out;


   offset += towerH*Sin(cellData.ThetaMin());

}// end RenderBarrelCell

//______________________________________________________________________________
void FWEveWriteToObjFile::RenderEndCapCell(const TEveCaloData::CellGeom_t &cellData, Float_t towerH, Float_t& offset, Float_t* pnts) const
{
   // Render an endcap cell.

   using namespace TMath;
   Float_t z1, r1In, r1Out, z2, r2In, r2Out;

   if (cellData.EtaMin()<0) z1  = endcapB_z;
   else  z1  = endcapF_z;

   z2    = z1 + TMath::Sign(towerH, cellData.EtaMin());

   r1In  = z1*Tan(cellData.ThetaMin());
   r2In  = z2*Tan(cellData.ThetaMin());
   r1Out = z1*Tan(cellData.ThetaMax());
   r2Out = z2*Tan(cellData.ThetaMax());

   Float_t cos2 = Cos(cellData.PhiMin());
   Float_t sin2 = Sin(cellData.PhiMin());
   Float_t cos1 = Cos(cellData.PhiMax());
   Float_t sin1 = Sin(cellData.PhiMax());

   //   Float_t box[24];
   //Float_t* pnts = box;
   // 0
   pnts[0] = r1In*cos1;
   pnts[1] = r1In*sin1;
   pnts[2] = z1;
   pnts += 3;
   // 1
   pnts[0] = r1In*cos2;
   pnts[1] = r1In*sin2;
   pnts[2] = z1;
   pnts += 3;
   // 2
   pnts[0] = r2In*cos2;
   pnts[1] = r2In*sin2;
   pnts[2] = z2;
   pnts += 3;
   // 3
   pnts[0] = r2In*cos1;
   pnts[1] = r2In*sin1;
   pnts[2] = z2;
   pnts += 3;
   //---------------------------------------------------
   // 4
   pnts[0] = r1Out*cos1;
   pnts[1] = r1Out*sin1;
   pnts[2] = z1;
   pnts += 3;
   // 5
   pnts[0] = r1Out*cos2;
   pnts[1] = r1Out*sin2;
   pnts[2] = z1;
   pnts += 3;
   // 6
   pnts[0] = r2Out*cos2;
   pnts[1] = r2Out*sin2;
   pnts[2] = z2;
   pnts += 3;
   // 7
   pnts[0] = r2Out*cos1;
   pnts[1] = r2Out*sin1;
   pnts[2] = z2;


   if (z1 > 0)
      offset += towerH * Cos(cellData.ThetaMin());
   else
      offset -= towerH * Cos(cellData.ThetaMin());

} // end RenderEndCapCell
