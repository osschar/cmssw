#ifndef Fireworks_Core_FWEveWriteToObjFile_h
#define Fireworks_Core_FWEveWriteToObjFile_h

#include "TEveCaloData.h"

#include <vector>

class TEveCalo3D;
class TGListTreeItem;
class TEveCalo3D;
class TEveTrack;


class FWEveWriteToObjFile
{
private:

  FWEveWriteToObjFile(const FWEveWriteToObjFile&);            // Not implemented
  FWEveWriteToObjFile& operator=(const FWEveWriteToObjFile&); // Not implemented

  void    CrossProduct(const Float_t a[3], const Float_t b[3], const Float_t c[3], Float_t out[3]) const;
  void    RenderGridEndCap() const;
  void    RenderGridBarrel() const;
  void    RenderGrid() const;
  void    RenderBarrelCell(const TEveCaloData::CellGeom_t &cell, Float_t towerH, Float_t& offset, Float_t* pnts) const;
  void    RenderEndCapCell(const TEveCaloData::CellGeom_t &cell, Float_t towerH, Float_t& offset, Float_t* pnts) const;

  TString give_track_vertecies(TEveTrack* akt_track);
  TString give_box_vetrecies_and_norm(Float_t* pnts);
  void write_tracks(std::ofstream &out_file, const TString& track_name, TGListTreeItem *akt_tree_item, int &n_track);
  void write_calo(std::ofstream &out_file, const TString& calo_name);
  void write_chambers(std::ofstream &out_file, const TString& chamber_name, TGListTreeItem *akt_tree_item, int &n_chamber);

  double scale;

  int n_vert;
  int n_norm;

protected:

  TEveCalo3D* fM;
  mutable std::vector<Float_t>     fOffset;

  TGListTreeItem *eventscene_3dtower;
  TGListTreeItem *tracks_from_3dtower;
  TGListTreeItem *ecal_from_3dtower;
  TGListTreeItem *muons_from_3dtower;
  Float_t transF_eta;
  Float_t transB_eta;
  Float_t endcapF_z;
  Float_t endcapB_z;

public:

  FWEveWriteToObjFile(TGListTreeItem*);

  void save_to_obj(const char* iName,
                   bool save_tracks    =true, bool save_calo        =true,
                   bool save_muon_track=true, bool save_muon_chamber=true);

  TGListTreeItem *my_FindSiblingByName(TGListTreeItem *item, const TString& name);
};


inline void FWEveWriteToObjFile::CrossProduct(const Float_t a[3], const Float_t b[3],
                                              const Float_t c[3], Float_t out[3]) const
{
   // Calculate cross-product.

   const Float_t v1[3] = { a[0] - c[0], a[1] - c[1], a[2] - c[2] };
   const Float_t v2[3] = { b[0] - c[0], b[1] - c[1], b[2] - c[2] };

   out[0] = v1[1] * v2[2] - v1[2] * v2[1];
   out[1] = v1[2] * v2[0] - v1[0] * v2[2];
   out[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

#endif
