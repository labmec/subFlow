//
//  TRMSpatialPropertiesMap.cpp
//  PZ
//
//  Created by omar duran on 21/06/2020.
//
//

#include "TRMSpatialPropertiesMap.h"
#include <limits>
#include <cassert>

TRMSpatialPropertiesMap::TRMSpatialPropertiesMap(){
    
    m_grid_coordinates.clear();
    m_properties.clear();
    
    m_kappa_default = 1.0e-7;
    m_phi_default = 0.1;
    
    m_n_SAMe_blocks.clear();
    m_size_SAMe_blocks.clear();
    m_SAMe_to_SPE.clear();
    m_n_spe_blocks.clear();
    m_size_spe_blocks.clear();
    m_spe_translation.clear();
    
    m_x_box.first = std::numeric_limits<REAL>::max();
    m_x_box.second = std::numeric_limits<REAL>::min();
    m_y_box.first = std::numeric_limits<REAL>::max();
    m_y_box.second = std::numeric_limits<REAL>::min();
    m_z_box.first = std::numeric_limits<REAL>::max();
    m_z_box.second = std::numeric_limits<REAL>::min();
    
}

TRMSpatialPropertiesMap::~TRMSpatialPropertiesMap(){
    
}

void TRMSpatialPropertiesMap::SampleKappaAndPhi(TPZManVector<REAL,3> &x, std::vector<REAL> &kappa_and_phi){
    if (m_map_type == ECartesianGrid) {
        SampleKappaAndPhiCartesianGrid(x,kappa_and_phi);
    }else if (m_map_type == ECornerPointGrid){
        SampleKappaAndPhiCornerGrid(x,kappa_and_phi);
    }else{
        assert(m_map_type == ENone);
    }
}



void TRMSpatialPropertiesMap::SetCartesianMeshData(std::vector<size_t> n_blocks, std::vector<REAL> size_blocks, std::string perm_data_name, std::string phi_data_name, std::vector<size_t> n_SAMe_blocks, std::vector<REAL> translation){
    
    m_map_type = ECartesianGrid;
    
    assert(n_blocks.size() == 3);
    
    assert(size_blocks.size() == 3);
    
    assert(n_SAMe_blocks.size() == 3);
    
    m_n_spe_blocks = n_blocks;
    m_size_spe_blocks = size_blocks;
    m_n_SAMe_blocks = n_SAMe_blocks;
    m_spe_translation = translation;
    
    std::ifstream stream_perm (perm_data_name.c_str());
    std::ifstream stream_phi (phi_data_name.c_str());
    
    size_t n_cells = m_n_spe_blocks[0] *  m_n_spe_blocks[1] * m_n_spe_blocks[2];
    m_properties.resize(n_cells);
    size_t c = 0;
    for (unsigned int k = 0; k < m_n_spe_blocks[2]; k++) {
        for (unsigned int j = 0; j < m_n_spe_blocks[1]; j++) {
            for (unsigned int i = 0; i < m_n_spe_blocks[0]; i++) {
                REAL k,phi;
                stream_perm >> k;
                stream_phi >> phi;
                m_properties[c] = {k,k,k,phi};
                c++;
            }
        }
    }
    
    assert(c == n_cells);
    
    BuildSAMeForCartesianMesh();
}

void TRMSpatialPropertiesMap::BuildSAMeForCartesianMesh(){
        
    auto substract = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] - b[0],a[1] - b[1],a[2] - b[2]};
    };
    
    auto add = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] + b[0],a[1] + b[1],a[2] + b[2]};
    };
    
    auto dot = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> REAL {
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    
    auto VoxelCenter = [] (unsigned int & i, unsigned int & j, unsigned int & k, REAL & dx, REAL & dy, REAL & dz) -> std::vector<REAL> {
            return {-dx/2. + dx*(i+1),-dy/2. + dy*(j+1),-dz/2. + dz*(k+1)};
    };
    
    auto IsVoxelMemberQ = [&substract,&dot] (REAL & dx, REAL & dy, REAL & dz, const std::vector<REAL> & center, const std::vector<REAL> & pt) -> bool {
        std::vector<REAL> e1 = {1,0,0};
        std::vector<REAL> e2 = {0,1,0};
        std::vector<REAL> e3 = {0,0,1};
        std::vector<REAL> d = substract(pt,center);
        bool checkxQ = std::fabs(dot(d,e1)) <= dx/2.0;
        bool checkyQ = std::fabs(dot(d,e2)) <= dy/2.0;
        bool checkzQ = std::fabs(dot(d,e3)) <= dz/2.0;
        return checkxQ && checkyQ && checkzQ;
    };
        
    REAL dx = m_size_spe_blocks[0];
    REAL dy = m_size_spe_blocks[1];
    REAL dz = m_size_spe_blocks[2];
    
    m_size_SAMe_blocks.resize(3);
    m_size_SAMe_blocks[0] = dx*m_n_spe_blocks[0]/m_n_SAMe_blocks[0];
    m_size_SAMe_blocks[1] = dy*m_n_spe_blocks[1]/m_n_SAMe_blocks[1];
    m_size_SAMe_blocks[2] = dz*m_n_spe_blocks[2]/m_n_SAMe_blocks[2];
    
    REAL sdx = m_size_SAMe_blocks[0];
    REAL sdy = m_size_SAMe_blocks[1];
    REAL sdz = m_size_SAMe_blocks[2];
    for (unsigned int sk = 0; sk < m_n_SAMe_blocks[2]; sk++) {
        for (unsigned int sj = 0; sj < m_n_SAMe_blocks[1]; sj++) {
            for (unsigned int si = 0; si < m_n_SAMe_blocks[0]; si++) {
                auto SAMe_voxcel_center = add(VoxelCenter(si,sj,sk,sdx,sdy,sdz),m_spe_translation);
                
                for (unsigned int k = 0; k < m_n_spe_blocks[2]; k++) {
                    for (unsigned int j = 0; j < m_n_spe_blocks[1]; j++) {
                        for (unsigned int i = 0; i < m_n_spe_blocks[0]; i++) {
                            auto spe_voxcel_center = add(VoxelCenter(i,j,k,dx,dy,dz),m_spe_translation);
                            bool isSAMeMemberQ = IsVoxelMemberQ(sdx,sdy,sdz,SAMe_voxcel_center,spe_voxcel_center);
                            if (isSAMeMemberQ) {
                                m_SAMe_to_SPE[{si,sj,sk}].push_back({i,j,k});
                            }
                        }
                    }
                }
            }
        }
    }
}

void TRMSpatialPropertiesMap::SampleKappaAndPhiCartesianGrid(TPZManVector<REAL,3> &x, std::vector<REAL> &kappa_and_phi){
    
    std::vector<REAL> pt = {x[0],x[1],x[2]};
    
    auto substract = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] - b[0],a[1] - b[1],a[2] - b[2]};
    };
    
    auto add = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] + b[0],a[1] + b[1],a[2] + b[2]};
    };
    
    auto dot = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> REAL {
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    
    auto VoxelCenter = [] (unsigned int & i, unsigned int & j, unsigned int & k, REAL & dx, REAL & dy, REAL & dz) -> std::vector<REAL> {
            return {-dx/2. + dx*(i+1),-dy/2. + dy*(j+1),-dz/2. + dz*(k+1)};
    };
    
    auto IsVoxelMemberQ = [&substract,&dot] (REAL & dx, REAL & dy, REAL & dz, const std::vector<REAL> & center, const std::vector<REAL> & pt) -> bool {
        std::vector<REAL> e1 = {1,0,0};
        std::vector<REAL> e2 = {0,1,0};
        std::vector<REAL> e3 = {0,0,1};
        std::vector<REAL> d = substract(pt,center);
        bool checkxQ = std::fabs(dot(d,e1)) <= dx/2.0;
        bool checkyQ = std::fabs(dot(d,e2)) <= dy/2.0;
        bool checkzQ = std::fabs(dot(d,e3)) <= dz/2.0;
        return checkxQ && checkyQ && checkzQ;
    };
        
    REAL dx = m_size_spe_blocks[0];
    REAL dy = m_size_spe_blocks[1];
    REAL dz = m_size_spe_blocks[2];
    
    REAL sdx = m_size_SAMe_blocks[0];
    REAL sdy = m_size_SAMe_blocks[1];
    REAL sdz = m_size_SAMe_blocks[2];
    
    unsigned int pos = -1;
    unsigned int i,j,k;
    for (unsigned int sk = 0; sk < m_n_SAMe_blocks[2]; sk++) {
        for (unsigned int sj = 0; sj < m_n_SAMe_blocks[1]; sj++) {
            for (unsigned int si = 0; si < m_n_SAMe_blocks[0]; si++) {
                auto SAMe_voxcel_center = add(VoxelCenter(si,sj,sk,sdx,sdy,sdz),m_spe_translation);
                
                bool isSAMeMemberQ = IsVoxelMemberQ(sdx,sdy,sdz,SAMe_voxcel_center,pt);
                
                if (isSAMeMemberQ) {
                    
                    std::map<std::vector<unsigned int>,std::vector<std::vector<unsigned int>>>::iterator voxels;
                    voxels = m_SAMe_to_SPE.find({si,sj,sk});
                    if (voxels == m_SAMe_to_SPE.end()){
                        continue;
                    }
                    for (auto voxcel : voxels->second) {
                        i = voxcel[0];
                        j = voxcel[1];
                        k = voxcel[2];
                        auto spe_voxcel_center = add(VoxelCenter(i,j,k,dx,dy,dz),m_spe_translation);
                        bool isSPEMemberQ = IsVoxelMemberQ(dx,dy,dz,spe_voxcel_center,pt);
                        if (isSPEMemberQ) {
                            pos = i + m_n_spe_blocks[0]*(j) + m_n_spe_blocks[0]*m_n_spe_blocks[1]*(k);
                            std::vector<REAL> chunk = m_properties[pos];
                            kappa_and_phi = chunk;
                            REAL mDTom2 = 0.986923e-15;
                            kappa_and_phi[0] *= mDTom2*1.0e6;
                            kappa_and_phi[1] *= mDTom2*1.0e6;
                            kappa_and_phi[2] *= mDTom2*1.0e6;
                            kappa_and_phi[3] += 1.0e-4; // Because some porosities are null
                            return;
                        }
                    }

                }
            }
        }
    }
    
    kappa_and_phi = {m_kappa_default,m_kappa_default,m_kappa_default,m_phi_default};
    return;
}

void TRMSpatialPropertiesMap::SampleKappaAndPhiCornerGrid(TPZManVector<REAL,3> &x, std::vector<REAL> &kappa_and_phi){
    
    std::vector<REAL> pt = {x[0],x[1],x[2]};
    
    auto substract = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] - b[0],a[1] - b[1],a[2] - b[2]};
    };
    
    auto cross = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {-(a[2]*b[1]) + a[1]*b[2],a[2]*b[0] - a[0]*b[2],-(a[1]*b[0]) + a[0]*b[1]};
    };
    
    auto dot = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> REAL {
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    
    auto sign = [] (const REAL & a) -> int {
            return (REAL(0.0) < a) - (a < REAL(0.0));
    };
    
    auto SameSide = [&substract,&cross,&dot,&sign] (const std::vector<REAL> & v1, const std::vector<REAL> & v2, const std::vector<REAL> & v3, const std::vector<REAL> & v4, const std::vector<REAL> & pt) -> bool {
        std::vector<REAL> normal = cross(substract(v2,v1),substract(v3,v1));
        REAL dot_v4 = dot(normal,substract(v4,v1));
        REAL dot_pt = dot(normal,substract(pt,v1));
        bool checkQ = sign(dot_v4) == sign(dot_pt);
        REAL tol = 1.0e-8;
        if (fabs(dot_v4) < tol) { // collapsed tet
            checkQ = false;
        }
        return checkQ;
    };
    
    auto IsTetrahedronMemberQ = [&SameSide] (const std::vector<REAL> & v1, const std::vector<REAL> & v2, const std::vector<REAL> & v3, const std::vector<REAL> & v4, const std::vector<REAL> & pt) -> bool {
        bool check1Q = SameSide(v1, v2, v3, v4, pt);
        bool check2Q = SameSide(v2, v3, v4, v1, pt);
        bool check3Q = SameSide(v3, v4, v1, v2, pt);
        bool check4Q = SameSide(v4, v1, v2, v3, pt);
        return check1Q && check2Q && check3Q && check4Q;
    };
    
    auto add = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] + b[0],a[1] + b[1],a[2] + b[2]};
    };
    
    auto VoxelCenter = [] (unsigned int & i, unsigned int & j, unsigned int & k, REAL & dx, REAL & dy, REAL & dz) -> std::vector<REAL> {
            return {-dx/2. + dx*(i+1),-dy/2. + dy*(j+1),-dz/2. + dz*(k+1)};
    };
    
    auto IsVoxelMemberQ = [&substract,&dot] (REAL & dx, REAL & dy, REAL & dz, const std::vector<REAL> & center, const std::vector<REAL> & pt) -> bool {
        std::vector<REAL> e1 = {1,0,0};
        std::vector<REAL> e2 = {0,1,0};
        std::vector<REAL> e3 = {0,0,1};
        std::vector<REAL> d = substract(pt,center);
        bool checkxQ = std::fabs(dot(d,e1)) <= dx/2.0;
        bool checkyQ = std::fabs(dot(d,e2)) <= dy/2.0;
        bool checkzQ = std::fabs(dot(d,e3)) <= dz/2.0;
        return checkxQ && checkyQ && checkzQ;
    };
    
    bool IsCellMemberQ = false;
    
    std::vector<REAL> translation(3);
    translation[0] = (m_x_box.first);
    translation[1] = (m_y_box.first);
    translation[2] = (m_z_box.first);
    
    REAL sdx = m_size_SAMe_blocks[0];
    REAL sdy = m_size_SAMe_blocks[1];
    REAL sdz = m_size_SAMe_blocks[2];
    for (unsigned int sk = 0; sk < m_n_SAMe_blocks[2]; sk++) {
        for (unsigned int sj = 0; sj < m_n_SAMe_blocks[1]; sj++) {
            for (unsigned int si = 0; si < m_n_SAMe_blocks[0]; si++) {
                auto SAMe_voxcel_center = add(VoxelCenter(si,sj,sk,sdx,sdy,sdz),translation);
                
                bool isSAMeMemberQ = IsVoxelMemberQ(sdx,sdy,sdz,SAMe_voxcel_center,pt);
                
                if (isSAMeMemberQ) {
                    
                    std::map<std::vector<unsigned int>,std::vector<unsigned int>>::iterator indexes;
                    indexes = m_SAMe_to_CornerGrid.find({si,sj,sk});
                    if (indexes == m_SAMe_to_CornerGrid.end()){
                        continue;
                    }
                    for (auto i : indexes->second) {
                        
                        std::vector<REAL> & cell_data = m_grid_coordinates[i];
                
                        std::vector<std::vector<REAL>> cell(8);
                        unsigned int c = 0;
                        for (unsigned int j = 0; j < 24; j += 3) {
                            cell[c] = {cell_data[j],cell_data[j+1],cell_data[j+2]};
                            c++;
                        }
                
                        std::vector<REAL> normal = cross(substract(cell[1],cell[0]),substract(cell[2],cell[0]));
                        bool isCollapsedQ = dot(normal,normal) < 1.0e-12;
                        if (isCollapsedQ) {
                            continue;
                        }
                
                        bool checktet1Q = IsTetrahedronMemberQ(cell[0],cell[1],cell[3],cell[4],pt);
                        bool checktet2Q = IsTetrahedronMemberQ(cell[3],cell[7],cell[4],cell[1],pt);
                        bool checktet3Q = IsTetrahedronMemberQ(cell[4],cell[5],cell[7],cell[1],pt);
                        bool checktet4Q = IsTetrahedronMemberQ(cell[1],cell[2],cell[3],cell[6],pt);
                        bool checktet5Q = IsTetrahedronMemberQ(cell[3],cell[6],cell[7],cell[1],pt);
                        bool checktet6Q = IsTetrahedronMemberQ(cell[5],cell[6],cell[7],cell[1],pt);
                
                        IsCellMemberQ = checktet1Q || checktet2Q || checktet3Q || checktet4Q || checktet5Q || checktet6Q;
                
                        if (IsCellMemberQ) {
                            std::vector<REAL> chunk = m_properties[i];
                            kappa_and_phi = chunk;
                            REAL mDTom2 = 0.986923e-15;
                            kappa_and_phi[0] *= mDTom2*1.0e6;
                            kappa_and_phi[1] *= mDTom2*1.0e6;
                            kappa_and_phi[2] *= mDTom2*1.0e6;
                            kappa_and_phi[3] += 1.0e-4; // Because some porosities are null
                            return;
                        }
                    }
                }
            }
        }
    }
    
    kappa_and_phi = {m_kappa_default,m_kappa_default,m_kappa_default,m_phi_default};
    return;
}

void TRMSpatialPropertiesMap::SetCornerGridMeshData(size_t n_cells, std::string corner_data_name, std::string props_data_name, std::vector<size_t> SAMe_blocks){
    
    m_map_type = ECornerPointGrid;
    
    assert(SAMe_blocks.size() == 3);
    
    m_n_SAMe_blocks = SAMe_blocks;
    
    std::ifstream stream_corners (corner_data_name.c_str());
    std::ifstream stream_props (props_data_name.c_str());
    REAL x,y,z,kx,ky,kz,phi;
    m_grid_coordinates.resize(n_cells);
    m_properties.resize(n_cells);
    
    for (unsigned int i = 0; i < n_cells; i++) {
        
        for (unsigned int j = 0; j < 24; j += 3) {
            stream_corners >> x;
            stream_corners >> y;
            stream_corners >> z;
            
            m_grid_coordinates[i].push_back(x);
            m_grid_coordinates[i].push_back(y);
            m_grid_coordinates[i].push_back(z);
            
            if (m_x_box.first > x) {
                m_x_box.first = x;
            }
            if (m_x_box.second < x) {
                m_x_box.second = x;
            }
            
            if (m_y_box.first > y) {
                m_y_box.first = y;
            }
            if (m_y_box.second < y) {
                m_y_box.second = y;
            }
            
            if (m_z_box.first > z) {
                m_z_box.first = z;
            }
            if (m_z_box.second < z) {
                m_z_box.second = z;
            }
            
        }
        
        stream_props >> kx;
        stream_props >> ky;
        stream_props >> kz;
        stream_props >> phi;
        
        m_properties[i].push_back(kx);
        m_properties[i].push_back(ky);
        m_properties[i].push_back(kz);
        m_properties[i].push_back(phi);
    }
    
    BuildSAMeForCornerGridMesh();
}



void TRMSpatialPropertiesMap::BuildSAMeForCornerGridMesh(){
    
    auto substract = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] - b[0],a[1] - b[1],a[2] - b[2]};
    };
    
    auto cross = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {-(a[2]*b[1]) + a[1]*b[2],a[2]*b[0] - a[0]*b[2],-(a[1]*b[0]) + a[0]*b[1]};
    };
    
    auto dot = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> REAL {
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    
    auto add = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] + b[0],a[1] + b[1],a[2] + b[2]};
    };
    
    auto VoxelCenter = [] (unsigned int & i, unsigned int & j, unsigned int & k, REAL & dx, REAL & dy, REAL & dz) -> std::vector<REAL> {
            return {-dx/2. + dx*(i+1),-dy/2. + dy*(j+1),-dz/2. + dz*(k+1)};
    };
    
    auto IsVoxelMemberQ = [&substract,&dot] (REAL & dx, REAL & dy, REAL & dz, const std::vector<REAL> & center, const std::vector<REAL> & pt) -> bool {
        std::vector<REAL> e1 = {1,0,0};
        std::vector<REAL> e2 = {0,1,0};
        std::vector<REAL> e3 = {0,0,1};
        std::vector<REAL> d = substract(pt,center);
        bool checkxQ = std::fabs(dot(d,e1)) <= dx/2.0;
        bool checkyQ = std::fabs(dot(d,e2)) <= dy/2.0;
        bool checkzQ = std::fabs(dot(d,e3)) <= dz/2.0;
        return checkxQ && checkyQ && checkzQ;
    };
    
    m_size_SAMe_blocks.resize(3);
    m_size_SAMe_blocks[0] = (m_x_box.second - m_x_box.first)/m_n_SAMe_blocks[0];
    m_size_SAMe_blocks[1] = (m_y_box.second - m_y_box.first)/m_n_SAMe_blocks[1];
    m_size_SAMe_blocks[2] = (m_z_box.second - m_z_box.first)/m_n_SAMe_blocks[2];
    
    std::vector<REAL> translation(3);
    translation[0] = (m_x_box.first);
    translation[1] = (m_y_box.first);
    translation[2] = (m_z_box.first);
    
    REAL sdx = m_size_SAMe_blocks[0];
    REAL sdy = m_size_SAMe_blocks[1];
    REAL sdz = m_size_SAMe_blocks[2];
    for (unsigned int sk = 0; sk < m_n_SAMe_blocks[2]; sk++) {
        for (unsigned int sj = 0; sj < m_n_SAMe_blocks[1]; sj++) {
            for (unsigned int si = 0; si < m_n_SAMe_blocks[0]; si++) {
                auto SAMe_voxcel_center = add(VoxelCenter(si,sj,sk,sdx,sdy,sdz),translation);

                for (unsigned int i = 0; i < m_grid_coordinates.size(); i++) {
                    std::vector<REAL> & cell_data = m_grid_coordinates[i];
                    
                    std::vector<std::vector<REAL>> cell(8);
                    unsigned int c = 0;
                    for (unsigned int j = 0; j < 24; j += 3) {
                        cell[c] = {cell_data[j],cell_data[j+1],cell_data[j+2]};
                        c++;
                    }
                    
                    std::vector<REAL> normal = cross(substract(cell[1],cell[0]),substract(cell[2],cell[0]));
                    bool isCollapsedQ = dot(normal,normal) < 1.0e-12;
                    if (isCollapsedQ) {
                        continue;
                    }
                    
                    std::vector<REAL> d1 = substract(cell[2],cell[0]);
                    std::vector<REAL> d2 = substract(cell[3],cell[1]);
                    REAL dh1 = std::sqrt(dot(d1,d1));
                    REAL dh2 = std::sqrt(dot(d2,d2));
                    std::vector<REAL> d1n = {d1[0]/dh1,d1[1]/dh1,d1[2]/dh1};
                    std::vector<REAL> d2n = {d2[0]/dh2,d2[1]/dh2,d2[2]/dh2};
                    std::vector<REAL> d1h = {d1n[0]*dh1/2.0,d1n[1]*dh1/2.0,d1n[2]*dh1/2.0};
                    std::vector<REAL> d2h = {d2n[0]*dh2/2.0,d2n[1]*dh2/2.0,d2n[2]*dh2/2.0};
                    std::vector<REAL> approx_center_d1 = add(d1h,cell[0]);
                    std::vector<REAL> approx_center_d2 = add(d2h,cell[1]);
                    
                    bool isSAMeMemberQ_1 = IsVoxelMemberQ(sdx,sdy,sdz,SAMe_voxcel_center,approx_center_d1);
                    bool isSAMeMemberQ_2 = IsVoxelMemberQ(sdx,sdy,sdz,SAMe_voxcel_center,approx_center_d2);
                    if (isSAMeMemberQ_1 || isSAMeMemberQ_2) {
                        m_SAMe_to_CornerGrid[{si,sj,sk}].push_back(i);
                    }
                    
                }
            }
        }
    }
}

void TRMSpatialPropertiesMap::Clear(){
    m_grid_coordinates.clear();
    m_properties.clear();
    m_n_SAMe_blocks.clear();
    m_size_SAMe_blocks.clear();
    m_SAMe_to_SPE.clear();
    m_n_spe_blocks.clear();
    m_size_spe_blocks.clear();
    m_spe_translation.clear();
    m_SAMe_to_CornerGrid.clear();
}
