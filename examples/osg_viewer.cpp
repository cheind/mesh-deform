/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#include "osg_viewer.h"

#include <osg/Group>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Camera>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osgFX/Scribe>
#include <osg/ShadeModel>
#include <osg/LightModel>
#include <osg/Material>

#include <unordered_set>

namespace deform {
    namespace example {
        
        struct OSGViewer::data {
            OSGViewer::Mesh *mesh;
            DeformCallback dc;
            osg::ref_ptr<osgViewer::Viewer> viewer;
            osg::ref_ptr<osg::Group> root;
            osg::ref_ptr<osg::Geode> geode;
        };
        
        class OSGUpdateCallback : public osg::NodeCallback
        {
        public:
            
            OSGUpdateCallback(OSGViewer::Mesh *mesh, OSGViewer::DeformCallback *dc)
            :_mesh(mesh), _dc(dc)
            {}
            
            virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
            {
                if ((*_dc)(*_mesh, nv->getFrameStamp()->getReferenceTime())) {
                    _mesh->update_normals();
                    
                    // Update vertex and normal buffer
                    osg::Geode *geode = static_cast<osg::Geode*>(node);
                    osg::Geometry *geo = static_cast<osg::Geometry*>(geode->getDrawable(0));
                    
                    osg::Vec3Array *vertices = static_cast<osg::Vec3Array*>(geo->getVertexArray());
                    osg::Vec3Array *normals = static_cast<osg::Vec3Array*>(geo->getNormalArray());
                    
                    for (auto v = _mesh->vertices_begin(); v != _mesh->vertices_end(); ++v) {
                        OSGViewer::Mesh::Point p = _mesh->point(*v);
                        vertices->at(v->idx()) = osg::Vec3(p[0], p[1], p[2]);
                        
                        OSGViewer::Mesh::Point n = _mesh->normal(*v);
                        normals->at(v->idx()) = osg::Vec3(n[0], n[1], n[2]);
                    }
                    
                    // Using VBOs, need to request re-upload.
                    vertices->dirty();
                    normals->dirty();
                }
                
                
                traverse(node,nv);
            }
            
        private:
            OSGViewer::Mesh *_mesh;
            OSGViewer::DeformCallback *_dc;
        };
        
        OSGViewer::OSGViewer(OSGViewer::Mesh &mesh)
        : OSGViewer(mesh, std::vector<int>(), std::vector<int>())
        {
        }
        
        OSGViewer::OSGViewer(Mesh &mesh, const std::vector<int> &anchors, const std::vector<int> &handles)
        :_data(new data())
        {
            _data->mesh = &mesh;
            _data->mesh->request_face_normals();
            _data->mesh->request_vertex_normals();
            _data->mesh->update_normals();
            
            setupScene(anchors, handles);
            
            _data->viewer = new osgViewer::Viewer();
            osg::DisplaySettings::instance()->setNumMultiSamples( 4 );
            _data->viewer->setUpViewInWindow(100, 100, 640, 360);
            _data->viewer->setSceneData(_data->root);
            _data->viewer->setThreadingModel(osgViewer::Viewer::SingleThreaded);
            _data->viewer->getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
            _data->viewer->getCamera()->setClearColor(osg::Vec4(1,1,1,1));
            
            _data->viewer->realize();
            
            typedef osgViewer::Viewer::Windows Windows;
            Windows windows;
            _data->viewer->getWindows(windows);
            for (auto w : windows) {
                w->setWindowName("Mesh-Deform https://github.com/cheind/mesh-deform");
            }
            
            _data->viewer->addEventHandler(new osgViewer::StatsHandler) ;
        }
        
        OSGViewer::~OSGViewer() {
            _data->mesh->release_face_normals();
            _data->mesh->release_vertex_normals();
        }
    
        
        void OSGViewer::setupScene(const std::vector<int> &anchors, const std::vector<int> &handles) {
            _data->root = new osg::Group();
            
            _data->geode = new osg::Geode();
            osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
            osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array();
            osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array();
        
            for (auto v = _data->mesh->vertices_begin(); v != _data->mesh->vertices_end(); ++v) {
                Mesh::Point p = _data->mesh->point(*v);
                vertices->push_back(osg::Vec3(p[0], p[1], p[2]));
                
                Mesh::Point n = _data->mesh->normal(*v);
                normals->push_back(osg::Vec3(n[0], n[1], n[2]));
            }
            
            
            osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array();
            colors->push_back(osg::Vec4(0.9, 0.9, 0.9, 1));
            colors->push_back(osg::Vec4(1.0, 0.0, 0.0, 1));
            colors->push_back(osg::Vec4(0.0, 1.0, 0.0, 1));
            
            osg::ref_ptr<osg::DrawElementsUInt> faces[3] = {
                new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0),
                new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0),
                new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0)
            };
            
            std::unordered_set<int> a_ids(anchors.begin(), anchors.end());
            std::unordered_set<int> h_ids(handles.begin(), handles.end());
            
            auto aend = a_ids.end();
            auto hend = h_ids.end();
            
            for (auto f = _data->mesh->faces_begin(); f != _data->mesh->faces_end(); ++f) {
                auto v = _data->mesh->fv_ccwbegin(*f);
                
                int id0 = v->idx(); ++v;
                int id1 = v->idx(); ++v;
                int id2 = v->idx();
                
                // Classify face
                bool isAnchor = a_ids.find(id0) != aend || a_ids.find(id1) != aend || a_ids.find(id2) != aend;
                bool isHandle = h_ids.find(id0) != hend || h_ids.find(id1) != hend || h_ids.find(id2) != hend;
                
                if (isAnchor) {
                    faces[1]->push_back(id0);
                    faces[1]->push_back(id1);
                    faces[1]->push_back(id2);
                } else if (isHandle) {
                    faces[2]->push_back(id0);
                    faces[2]->push_back(id1);
                    faces[2]->push_back(id2);
                } else {
                    faces[0]->push_back(id0);
                    faces[0]->push_back(id1);
                    faces[0]->push_back(id2);
                }
            }
            
            geometry->setVertexArray(vertices);
            geometry->setNormalArray(normals, osg::Array::BIND_PER_VERTEX);
            geometry->addPrimitiveSet(faces[0]);
            geometry->addPrimitiveSet(faces[1]);
            geometry->addPrimitiveSet(faces[2]);
            geometry->setUseVertexBufferObjects(true);
            geometry->setUseDisplayList(false);
            geometry->setColorArray(colors);
            geometry->setColorBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);
            
            _data->geode->addDrawable(geometry);
            _data->geode->setDataVariance(osg::Object::DYNAMIC);

            osg::StateSet* state = _data->geode->getOrCreateStateSet();
            osg::ref_ptr<osg::ShadeModel> sm = new osg::ShadeModel();
            sm->setMode(osg::ShadeModel::SMOOTH);
            state->setAttributeAndModes(sm, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
            
            osg::ref_ptr<osg::LightModel> lm = new osg::LightModel();
            lm->setTwoSided(true);
            state->setAttributeAndModes(lm, osg::StateAttribute::ON );
            
            osg::ref_ptr<osg::Material> mat = new osg::Material;
            mat->setColorMode(osg::Material::ColorMode::DIFFUSE);
            mat->setSpecular(osg::Material::FRONT, osg::Vec4f(0.8f, 0.8f, 0.8f, 1.0f));
            mat->setShininess(osg::Material::FRONT, 6.0f);
            state->setAttributeAndModes(mat, osg::StateAttribute::ON );
            
            state->setMode(GL_NORMALIZE, osg::StateAttribute::ON);
            
            _data->geode->setStateSet(state);
                        
            osg::ref_ptr<osgFX::Scribe> scribe = new osgFX::Scribe();
            scribe->setWireframeColor(osg::Vec4d(0.2,0.2,0.2,1));
            scribe->setWireframeLineWidth(0.5f);
            scribe->addChild(_data->geode);
            
            _data->root->addChild(scribe);
        }
        
        void OSGViewer::onFrame(DeformCallback dc) {
            _data->dc = dc;
            _data->geode->setUpdateCallback(new OSGUpdateCallback(_data->mesh, &_data->dc));
        }
        
        void OSGViewer::run() {
            _data->viewer->run();
        }
        
    }

}