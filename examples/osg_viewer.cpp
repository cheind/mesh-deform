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

namespace deform {
    namespace example {
        
        struct OSGViewer::data {
            OSGViewer::Mesh *mesh;
            DeformCallback dc;
            osg::ref_ptr<osgViewer::Viewer> viewer;
            osg::ref_ptr<osg::Group> root;
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
        
        OSGViewer::OSGViewer(int argc, char **argv, OSGViewer::Mesh &mesh, DeformCallback dc)
            : _data(new data())
        {
            _data->mesh = &mesh;
            _data->mesh->request_face_normals();
            _data->mesh->request_vertex_normals();
            _data->mesh->update_normals();
            _data->dc = dc;
            
            setupScene();
            
            osg::ArgumentParser arguments(&argc, argv);
            _data->viewer = new osgViewer::Viewer(arguments);
            _data->viewer->setUpViewInWindow(100, 100, 640, 360);
            //osg::DisplaySettings::instance()->setNumMultiSamples(8);
            _data->viewer->setSceneData(_data->root);
            _data->viewer->setThreadingModel(osgViewer::Viewer::SingleThreaded);
            _data->viewer->getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
            
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
        
        void OSGViewer::setupScene() {
            _data->root = new osg::Group();
            
            osg::ref_ptr<osg::Geode> geode = new osg::Geode();
            osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry();
            osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array();
            osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array();
        
            for (auto v = _data->mesh->vertices_begin(); v != _data->mesh->vertices_end(); ++v) {
                Mesh::Point p = _data->mesh->point(*v);
                vertices->push_back(osg::Vec3(p[0], p[1], p[2]));
                
                Mesh::Point n = _data->mesh->normal(*v);
                normals->push_back(osg::Vec3(n[0], n[1], n[2]));
            }
            
            osg::ref_ptr<osg::DrawElementsUInt> faces = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
            for (auto f = _data->mesh->faces_begin(); f != _data->mesh->faces_end(); ++f) {
                auto v = _data->mesh->fv_ccwbegin(*f);
                
                faces->push_back(v->idx()); ++v;
                faces->push_back(v->idx()); ++v;
                faces->push_back(v->idx());
                
                Mesh::Point p = _data->mesh->normal(*f);
                normals->push_back(osg::Vec3(p[0], p[1], p[2]));
            }
            
            geometry->setVertexArray(vertices);
            geometry->setNormalArray(normals, osg::Array::BIND_PER_VERTEX);
            geometry->addPrimitiveSet(faces);
            geometry->setUseVertexBufferObjects(true);
            geometry->setUseDisplayList(false);
            
            geode->addDrawable(geometry);
            geode->setDataVariance(osg::Object::DYNAMIC);
            geode->setUpdateCallback(new OSGUpdateCallback(_data->mesh, &_data->dc));

            osg::StateSet* state = geode->getOrCreateStateSet();
            osg::ShadeModel* sm = new osg::ShadeModel();
            sm->setMode(osg::ShadeModel::FLAT);
            state->setAttributeAndModes(sm, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
            geode->setStateSet(state);
                        
            osg::ref_ptr<osgFX::Scribe> scribe = new osgFX::Scribe();
            scribe->setWireframeColor(osg::Vec4d(0,0,0,1));
            scribe->setWireframeLineWidth(1.f);
            scribe->addChild(geode);
            
            _data->root->addChild(scribe);
        }
        
        void OSGViewer::run() {
            _data->viewer->run();
        }
        
    }

}