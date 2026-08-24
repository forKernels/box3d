// SPDX-FileCopyrightText: 2026 The Fantastic Planet
// SPDX-License-Identifier: MIT
//
// Coupling API: batch body snapshots, batch impulses, shape geometry, signed closest points.
// Written for the forKernels fork; intended to be upstreamable.

#include "test_macros.h"

#include "box3d/box3d.h"
#include "box3d/math_functions.h"
#include "box3d/types.h"

#include <math.h>

static int BodySnapshotsTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, 0.0f, 0.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );

	// dynamic unit sphere, density 1 -> mass 4/3 pi, moved and spun so every field is nontrivial
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = (b3Pos){ 1.0, 2.0, 3.0 };
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = 1.0f;
	b3Sphere sphere = { { 0.0f, 0.0f, 0.0f }, 1.0f };
	b3CreateSphereShape( bodyId, &shapeDef, &sphere );
	b3Body_SetLinearVelocity( bodyId, (b3Vec3){ 0.5f, -0.25f, 0.125f } );
	b3Body_SetAngularVelocity( bodyId, (b3Vec3){ 0.0f, 0.0f, 2.0f } );

	// a static body too, and one invalid id
	b3BodyDef staticDef = b3DefaultBodyDef();
	staticDef.type = b3_staticBody;
	b3BodyId staticId = b3CreateBody( worldId, &staticDef );
	b3BodyId ids[3] = { bodyId, staticId, b3_nullBodyId };
	b3BodySnapshot snap[3];
	int n = b3Body_GetSnapshots( ids, 3, snap );
	ENSURE( n == 3 );

	// dynamic body: every field equals the scalar accessor
	ENSURE( snap[0].isValid == true );
	ENSURE( snap[0].type == b3_dynamicBody );
	ENSURE_SMALL( snap[0].origin.x - 1.0, 1e-12 );
	ENSURE_SMALL( snap[0].origin.y - 2.0, 1e-12 );
	ENSURE_SMALL( snap[0].origin.z - 3.0, 1e-12 );
	ENSURE_SMALL( snap[0].rotation.s - 1.0f, 1e-6f );
	ENSURE_SMALL( snap[0].linearVelocity.x - 0.5f, 1e-6f );
	ENSURE_SMALL( snap[0].linearVelocity.y + 0.25f, 1e-6f );
	ENSURE_SMALL( snap[0].angularVelocity.z - 2.0f, 1e-6f );
	ENSURE_SMALL( snap[0].invMass - b3Body_GetInverseMass( bodyId ), 1e-7f );
	b3Matrix3 I = b3Body_GetLocalRotationalInertia( bodyId );
	ENSURE_SMALL( snap[0].localInertia.cx.x - I.cx.x, 1e-6f );
	ENSURE_SMALL( snap[0].localInertia.cy.y - I.cy.y, 1e-6f );
	ENSURE_SMALL( snap[0].localInertia.cz.z - I.cz.z, 1e-6f );
	// sphere: I = 2/5 m r^2, so I^-1 diag = 1/I
	ENSURE_SMALL( snap[0].localInvInertia.cx.x * I.cx.x - 1.0f, 1e-5f );
	ENSURE_SMALL( snap[0].localInvInertia.cy.x, 1e-6f );
	ENSURE_SMALL( snap[0].localCenter.x, 1e-6f );

	// static body: valid, zero inverse mass and inertia
	ENSURE( snap[1].isValid == true );
	ENSURE( snap[1].type == b3_staticBody );
	ENSURE( snap[1].invMass == 0.0f );
	ENSURE( snap[1].localInvInertia.cx.x == 0.0f && snap[1].localInvInertia.cy.y == 0.0f && snap[1].localInvInertia.cz.z == 0.0f );

	// invalid id: isValid false, zeroed
	ENSURE( snap[2].isValid == false );
	ENSURE( snap[2].invMass == 0.0f );

	b3DestroyWorld( worldId );
	return 0;
}

static int BodyImpulsesTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, 0.0f, 0.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = 1.0f;
	b3Sphere sphere = { { 0.0f, 0.0f, 0.0f }, 1.0f };
	b3CreateSphereShape( bodyId, &shapeDef, &sphere );
	float mass = b3Body_GetMass( bodyId );
	b3Matrix3 I = b3Body_GetLocalRotationalInertia( bodyId );

	// linear impulse P -> dv = P/m ; angular impulse L about z -> dw = L / Izz
	b3Vec3 lin[1] = { { 2.0f * mass, 0.0f, 0.0f } };
	b3Vec3 ang[1] = { { 0.0f, 0.0f, 3.0f * I.cz.z } };
	b3Body_ApplyImpulses( &bodyId, 1, lin, ang, true );
	b3Vec3 v = b3Body_GetLinearVelocity( bodyId );
	b3Vec3 w = b3Body_GetAngularVelocity( bodyId );
	ENSURE_SMALL( v.x - 2.0f, 1e-5f );
	ENSURE_SMALL( v.y, 1e-6f );
	ENSURE_SMALL( w.z - 3.0f, 1e-4f );

	// null arrays are zero impulses: nothing changes
	b3Body_ApplyImpulses( &bodyId, 1, NULL, NULL, true );
	b3Vec3 v2 = b3Body_GetLinearVelocity( bodyId );
	ENSURE_SMALL( v2.x - v.x, 1e-7f );

	// a static body is skipped, an invalid id is skipped, no crash
	b3BodyDef staticDef = b3DefaultBodyDef();
	staticDef.type = b3_staticBody;
	b3BodyId staticId = b3CreateBody( worldId, &staticDef );
	b3BodyId ids[2] = { staticId, b3_nullBodyId };
	b3Vec3 lin2[2] = { { 1.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f } };
	b3Body_ApplyImpulses( ids, 2, lin2, NULL, true );
	b3Vec3 vs = b3Body_GetLinearVelocity( staticId );
	ENSURE( vs.x == 0.0f );

	b3DestroyWorld( worldId );
	return 0;
}

int CouplingTest( void )
{
	RUN_SUBTEST( BodySnapshotsTest );
	RUN_SUBTEST( BodyImpulsesTest );
	return 0;
}
