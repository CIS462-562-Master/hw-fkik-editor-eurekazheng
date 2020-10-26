#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	AJoint *rootNode = m_pSkeleton->getRootNode();
	vec3 rootPos = m_Guide.getGlobalTranslation() + m_Guide.getGlobalRotation() * rootNode->getGlobalTranslation();

	// 2.	Set the y component of the guide position to 0
	rootPos[1] = 0;
	m_Guide.setGlobalTranslation(rootPos);

	// 3.	Set the global rotation of the guide joint towards the guideTarget
	guideTargetPos[1] = 0.f;
	vec3 origForward = vec3(0.f, 0.f, 1.f);
	vec3 newForward = guideTargetPos - m_Guide.getGlobalTranslation();
	vec3 axis = origForward ^ newForward;
	double angle = acos((origForward * newForward) / (origForward.Length() * newForward.Length()));
	mat3 rotM;
	rotM.FromAxisAngle(axis, angle);
	m_Guide.setGlobalRotation(rotM);
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space

	// 1.	Update the local translation of the root based on the left height and the right height
	AJoint *root = m_pSkeleton->getRootNode();
	vec3 rootPos = root->getGlobalTranslation();
	rootPos[1] += (leftHeight + rightHeight) / 2;
	root->setLocalTranslation(rootPos);

	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK 
	vec3 origUp = vec3(0.f, 1.f, 0.f);

	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal
		ATarget target;
		vec3 leftFootPos = leftFoot->getGlobalTranslation();
		leftFootPos[1] = leftHeight;
		target.setGlobalTranslation(leftFootPos);
		m_IKController->IKSolver_Limb(leftFoot->getID(), target);

		vec3 leftLocalNormal = leftFoot->getLocal2Global().m_rotation.Inverse() * leftNormal;
		vec3 axis = origUp ^ leftLocalNormal;
		double angle = acos((origUp * leftLocalNormal) / (origUp.Length() * leftLocalNormal.Length()));
		mat3 rotM;
		rotM.FromAxisAngle(axis, angle);
		leftFoot->setLocalRotation(rotM);
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal
		ATarget target;
		vec3 rightFootPos = rightFoot->getGlobalTranslation();
		rightFootPos[1] = rightHeight;
		target.setGlobalTranslation(rightFootPos);
		m_IKController->IKSolver_Limb(rightFoot->getID(), target);

		vec3 rightLocalNormal = rightFoot->getLocal2Global().m_rotation.Inverse() * rightNormal;
		vec3 axis = origUp ^ rightLocalNormal;
		double angle = acos((origUp * rightLocalNormal) / (origUp.Length() * rightLocalNormal.Length()));
		mat3 rotM;
		rotM.FromAxisAngle(axis, angle);
		rightFoot->setLocalRotation(rotM);
	}
	m_pSkeleton->update();
}
