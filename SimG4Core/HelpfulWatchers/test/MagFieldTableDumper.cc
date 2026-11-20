// MagFieldTableDumper: samples the CMSSW magnetic field (IdealMagneticFieldRecord, the field that the Geant4
// simulation wraps) on cylindrical grids and writes it as a flat ROOT tree, for REve's field model.
//
// Tree "field", one entry per grid point: grid (index into the "grids" VPSet), r, phi, z (cm, rad) and br, bphi,
// bz (tesla), all double. Each grid is r in [rMin, rMax] step dr, z in [zMin, zMax] step dz, at each phi of phiDeg.
// Also written, as TNamed: "config" (the "label" parameter, for the global tag and field configuration),
// "nominalValue" (MagneticField::nominalValue(), kGauss), "grids" (the grid definitions as text).

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TFile.h"
#include "TNamed.h"
#include "TTree.h"

class MagFieldTableDumper : public edm::one::EDAnalyzer<> {
public:
  explicit MagFieldTableDumper(edm::ParameterSet const& config)
      : fieldToken_(esConsumes()),
        fileName_(config.getParameter<std::string>("fileName")),
        label_(config.getParameter<std::string>("label")),
        grids_(config.getParameter<std::vector<edm::ParameterSet>>("grids")) {}

  void analyze(edm::Event const&, edm::EventSetup const& setup) override {
    if (done_)
      return;
    done_ = true;
    auto const& field = setup.getData(fieldToken_);

    TFile file(fileName_.c_str(), "RECREATE");
    TTree tree("field", "magnetic field on cylindrical grids: r, phi, z [cm, rad]; br, bphi, bz [T]");
    int grid = 0;
    double r = 0, phi = 0, z = 0, br = 0, bphi = 0, bz = 0;
    tree.Branch("grid", &grid);
    tree.Branch("r", &r);
    tree.Branch("phi", &phi);
    tree.Branch("z", &z);
    tree.Branch("br", &br);
    tree.Branch("bphi", &bphi);
    tree.Branch("bz", &bz);

    std::ostringstream desc;
    for (size_t g = 0; g < grids_.size(); ++g) {
      auto const& p = grids_[g];
      const double rMin = p.getParameter<double>("rMin"), rMax = p.getParameter<double>("rMax");
      const double dr = p.getParameter<double>("dr");
      const double zMin = p.getParameter<double>("zMin"), zMax = p.getParameter<double>("zMax");
      const double dz = p.getParameter<double>("dz");
      const auto phis = p.getParameter<std::vector<double>>("phiDeg");
      const int nr = std::lround((rMax - rMin) / dr), nz = std::lround((zMax - zMin) / dz);
      desc << "grid " << g << ": r " << rMin << ".." << rMax << " step " << dr << ", z " << zMin << ".." << zMax
           << " step " << dz << " (cm), phi [deg]";
      for (double ph : phis)
        desc << " " << ph;
      desc << "\n";
      grid = g;
      for (double phDeg : phis) {
        phi = phDeg * M_PI / 180.;
        const double c = std::cos(phi), s = std::sin(phi);
        for (int ir = 0; ir <= nr; ++ir) {
          r = rMin + ir * dr;
          for (int iz = 0; iz <= nz; ++iz) {
            z = zMin + iz * dz;
            const GlobalVector b = field.inTesla(GlobalPoint(r * c, r * s, z));
            br = b.x() * c + b.y() * s;
            bphi = -b.x() * s + b.y() * c;
            bz = b.z();
            tree.Fill();
          }
        }
      }
    }
    tree.Write();
    TNamed("config", label_.c_str()).Write();
    TNamed("nominalValue", std::to_string(field.nominalValue()).c_str()).Write();
    TNamed("grids", desc.str().c_str()).Write();
    file.Close();
  }

private:
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;
  const std::string fileName_, label_;
  const std::vector<edm::ParameterSet> grids_;
  bool done_ = false;
};

DEFINE_FWK_MODULE(MagFieldTableDumper);
