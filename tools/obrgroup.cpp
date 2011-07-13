#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/alias.h>
#include <openbabel/mcdlutil.h>

using namespace OpenBabel;

struct RGroup
{
  RGroup() : mol(0), atom(0), doubleBond(false) {}

  OBMol *mol; // the R group molecule
  OBAtom *atom; // the atom to attach to
  std::string rrr; // e.g. R2
  bool doubleBond; // is the bond to R# a double bond
};

void EnumerateProducts(OBMol *scaffold, std::map<std::string, OBAtom*> &scaffoldRMap, 
    std::map<std::string, std::vector<RGroup> > &rGroupMap,
    std::map<std::string, std::vector<RGroup> >::iterator iter, OBConversion &conv)
{
  std::string rrr = iter->first;
  std::vector<RGroup> &rGroups = iter->second;
  iter++;
  for (std::size_t i = 0; i < rGroups.size(); ++i) {
    OBMol scaffoldCopy(*scaffold);
    if (scaffoldRMap.find(rrr) == scaffoldRMap.end()) {
      std::cout << "Warning: Could not find substituent R# in scaffold. Ignoring substituent." << std::endl;
      break;
    }
    OBAtom *scaffoldAtom = scaffoldCopy.GetAtom(scaffoldRMap[rrr]->GetIdx());

    scaffoldCopy += *rGroups[i].mol;
    OBAtom *rGroupAtom = scaffoldCopy.GetAtom(scaffold->NumAtoms() + rGroups[i].atom->GetIdx());

    scaffoldCopy.AddBond(scaffoldAtom->GetIdx(), rGroupAtom->GetIdx(), 1);

    if (iter != rGroupMap.end())
      EnumerateProducts(&scaffoldCopy, scaffoldRMap, rGroupMap, iter, conv);
    else {
      generateDiagram(&scaffoldCopy);
      conv.Write(&scaffoldCopy);
      std::cout << "product..." << std::endl;
    }
  }
}

int main(int argc, char **argv)
{
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <input> <output>" << std::endl;
    return 0;
  }

  std::string input = argv[1];
  std::string output = argv[2];

  std::ifstream ifs(input.c_str());
  std::ofstream ofs(output.c_str());

  OBConversion conv(&ifs, &ofs);
  OBFormat *inFormat = conv.FormatFromExt(input);
  OBFormat *outFormat = conv.FormatFromExt(output);

  if (!inFormat || !conv.SetInFormat(inFormat)) {
    std::cout << "Could not find the input file format." << std::endl;
    return 0;
  }

  if (!outFormat || !conv.SetOutFormat(outFormat)) {
    std::cout << "Could not find the output file format." << std::endl;
    return 0;
  }

  OBMol scaffold;

  if (!conv.Read(&scaffold)) {
    std::cout << "Could not read the scaffold molecule." << std::endl;
    return 0;
  }

  // map R# to neighbor atom
  std::map<std::string, OBAtom*> scaffoldRMap;

  std::vector<OBAtom*> atomsToDelete;
  FOR_ATOMS_OF_MOL (atom, scaffold) {
    if (atom->HasData(AliasDataType)) {
      OBAtom *nbrAtom = 0;
      FOR_NBORS_OF_ATOM (nbr, &*atom) {
        nbrAtom = &*nbr;
      }
      AliasData *ad = static_cast<AliasData*>(atom->GetData(AliasDataType));
      scaffoldRMap[ad->GetAlias()] = nbrAtom;
      std::cout << ad->GetAlias() << " -> " << nbrAtom->GetIdx() << std::endl;
      atomsToDelete.push_back(&*atom);
    }
  }

  // Delete the R# dummy atoms
  for (std::size_t i = 0; i < atomsToDelete.size(); ++i)
    scaffold.DeleteAtom(atomsToDelete[i]);


  std::cout << "Searching for R# substituents..." << std::endl;

  // Find the R# groups
  std::map<std::string, std::vector<RGroup> > rGroupMap;
  while (true) {
    OBMol *mol = new OBMol;
    if (!conv.Read(mol)) {
      delete mol;
      break;
    }

    OBAtom *atomToDelete = 0;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (atom->HasData(AliasDataType)) {
        if (atomToDelete) {
          std::cout << "Error: substituent has more than one R#" << std::endl;
          break;
        }
        OBAtom *nbrAtom = 0;
        FOR_NBORS_OF_ATOM (nbr, &*atom) {
          nbrAtom = &*nbr;
        }
        AliasData *ad = static_cast<AliasData*>(atom->GetData(AliasDataType));

        if (rGroupMap.find(ad->GetAlias()) == rGroupMap.end())
          rGroupMap[ad->GetAlias()] = std::vector<RGroup>();
        std::vector<RGroup> &rGroups = rGroupMap[ad->GetAlias()];

        rGroups.push_back(RGroup());
        RGroup &rGroup = rGroups.back();

        rGroup.mol = mol;
        rGroup.rrr = ad->GetAlias();
        rGroup.atom = nbrAtom;

        atomToDelete = &*atom;
        std::cout << ad->GetAlias() << " -> " << nbrAtom->GetIdx() << std::endl;
      }
    }

    if (!atomToDelete) {
      delete mol;
      continue;
    }

    mol->DeleteAtom(atomToDelete);
  }

  std::cout << "Enumerating products..." << std::endl;

  EnumerateProducts(&scaffold, scaffoldRMap, rGroupMap, rGroupMap.begin(), conv);


}
