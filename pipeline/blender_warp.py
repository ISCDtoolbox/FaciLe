import os
import bpy
import sys

def shrinkwrap(obj, source=None, iterations=4):
    if source is None:
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=4, location=(0.5, 0.5, 0.5))
        source = bpy.context.active_object
    for i in range(iterations):
        print ("Iteration " + str(i))
        bpy.context.scene.objects.active = source
        bpy.ops.object.modifier_add(type='SHRINKWRAP')
        bpy.context.object.modifiers["Shrinkwrap"].target = obj
        bpy.context.object.modifiers["Shrinkwrap"].use_keep_above_surface = True
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Shrinkwrap")
        bpy.ops.object.modifier_add(type='REMESH')
        bpy.context.object.modifiers["Remesh"].octree_depth = 8
        bpy.context.object.modifiers["Remesh"].mode = 'SMOOTH'
        bpy.context.object.modifiers["Remesh"].use_smooth_shade = True
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Remesh")
    return source

if __name__=="__main__":
    """
    argv[1] = cage   file
    argv[2] = target file
    argv[3] = output file
    """
    argv = sys.argv[sys.argv.index("--") + 1:]
    bpy.ops.import_mesh.mesh('INVOKE_DEFAULT', filepath=argv[0])
    cage   = bpy.context.active_object
    bpy.ops.import_mesh.mesh('INVOKE_DEFAULT', filepath=argv[1])
    target = bpy.context.active_object

    warped = shrinkwrap(target, source=cage)
    bpy.ops.object.select_all(action='DESELECT')
    warped.select=True
    bpy.context.scene.objects.active = warped

    for ob in bpy.data.objects:
        for vgroup in ob.vertex_groups:
            ob.vertex_groups.remove(vgroup)

    bpy.ops.export_mesh.mesh('INVOKE_DEFAULT', filepath=argv[2])
