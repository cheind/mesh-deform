
# Note, mesh needs to be in object mode
verts = [i.index for i in bpy.context.active_object.data.vertices if i.select]
print(','.join(str(x) for x in verts))
