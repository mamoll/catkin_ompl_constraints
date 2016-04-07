/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2014, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Ryan Luna */

#include <algorithm> // for std::swap
#include "ompl/base/ConstrainedSpaceInformation.h"
#include "ompl/geometric/PathGeometric.h"

ompl::base::ConstrainedSpaceInformation::ConstrainedSpaceInformation(const StateSpacePtr &space)
    : ompl::base::SpaceInformation(space)
{
}

void ompl::base::ConstrainedSpaceInformation::setConstraintInformation(const ConstraintInformationPtr& ci)
{
    ci_ = ci;
}

const ompl::base::ConstraintInformationPtr& ompl::base::ConstrainedSpaceInformation::getConstraintInformation() const
{
    return ci_;
}

bool ompl::base::ConstrainedSpaceInformation::constrainedExtend(const base::State *a, const base::State *b,
                                                 std::vector<base::State*>& result) const
{
    // A semi-faithful implementation from the paper
    // Instead of using vector operations, this implementation
    // will accomplish a similar extension using interpolation in
    // the state space

    // Assuming a and b are both valid and on constraint manifold
    base::StateSpacePtr ss = getStateSpace();
    const base::State* previous = a;

    // number of discrete steps between a and b in the state space
    int n = ss->validSegmentCount(a, b);

    if (n == 0) // don't divide by zero
        return true;

    double dist = distance(a, b);
    double delta = dist / n; // This is the step size that we will take during extension

    while (true)
    {
        // The distance to travel is less than our step size.  Just declare victory
        if (dist < (delta + std::numeric_limits<double>::epsilon()))
        {
            result.push_back(cloneState(b));
            return true;
        }

        // Compute the parametrization for interpolation
        double t = delta / dist;

        base::State* scratchState = allocState();
        ss->interpolate(previous, b, t, scratchState);

        // Project new state onto constraint manifold.  Make sure the new state is valid
        // and that it has not deviated too far from where we started
        if (!ci_->project(scratchState) ||
            !isValid(scratchState) ||
            distance(previous, scratchState) > 2.0*delta)
        {
            freeState(scratchState);
            break;
        }

        // Check for divergence.  Divergence is declared if we are no closer to b
        // than before projection
        double newDist = distance(scratchState, b);
        if (newDist >= dist)
        {
            // Since we already collision checked this state, we might as well keep it
            result.push_back(scratchState);
            break;
        }
        dist = newDist;

        // No divergence; getting closer.  Store the new state
        result.push_back(scratchState);
        previous = scratchState;
    }
    return false;
}
bool ompl::base::ConstrainedSpaceInformation::shortcutPath(ompl::geometric::PathGeometric *path,
	const base::PlannerTerminationCondition &ptc) const
{
    std::vector<base::State*>& states = path->getStates();

    unsigned int count = 0;
    unsigned int maxCount = 10;

    while(!ptc && count < maxCount && states.size() > 3)
    {
        int i = rng_.uniformInt(0, states.size()-2);
        int j;
        do
        {
            j = rng_.uniformInt(0, states.size()-1);
        } while (abs(i-j) < 2); // make sure the difference between i and j is at least two

        if (i > j)
            std::swap(i, j);

        base::State* a = states[i];
        base::State* b = states[j];
        std::vector<base::State*> shortcut;

        bool foundShortcut = false;
        if (constrainedExtend(a, b, shortcut) && shortcut.size() > 1)
        {
            // see if shortcut is actually shorter
            double shortcutDist = 0.0;
            for(size_t k = 1; k < shortcut.size(); ++k)
                shortcutDist += distance(shortcut[k-1], shortcut[k]);

            double pathDist = 0.0;
            for(int k = i+1; k < j; ++k)
                pathDist += distance(states[k-1], states[k]);

            if (shortcutDist < pathDist)
            {
                // Delete states between [i+1, j]
                for(int k = i+1; k < j+1; ++k)
                    freeState(states[k]);
                states.erase(states.begin() + i+1, states.begin() + j+1);

                // Inserting new shortcut
                states.insert(states.begin() + i+1, shortcut.begin(), shortcut.end());
                foundShortcut = true;
            }
        }

        if (!foundShortcut)
        {
            for(size_t k = 0; k < shortcut.size(); ++k)
                freeState(shortcut[k]);
            count++;
        }
        else
            count = 0;
    }

    return true;
}

bool ompl::base::ConstrainedSpaceInformation::projectPath(
	const geometric::PathGeometric &inpath, geometric::PathGeometric &outpath,
	bool waypointsValid) const
{
	unsigned int numStates = inpath.getStateCount();
	bool result = true;
	std::vector<State*> projectedWaypoints;

    if (waypointsValid)
        for (unsigned int i = 0; i < inpath.getStateCount(); ++i)
            assert(isValid(inpath.getState(i)));

	if (!waypointsValid)
	{
        OMPL_WARN("Projecting a path with %d waypoints onto constraint manifold", numStates);
		base::State* scratchState = allocState();
		for (unsigned int i = 0; i < numStates; ++i)
		{
            bool p, v;
			copyState(scratchState, inpath.getState(i));
			if (!(p=ci_->project(scratchState)) || !(v=isValid(scratchState)))
			{
                OMPL_WARN("Projection of waypoint %d failed! projection=%d, validity=%d",
                    i, p, v);
				freeState(scratchState);
				for (unsigned int j = 0; j < projectedWaypoints.size(); ++j)
					freeState(projectedWaypoints[j]);
				return false;
			}
			projectedWaypoints.push_back(cloneState(scratchState));
		}
		freeState(scratchState);
	}

	const std::vector<State*>& instates = waypointsValid
		? const_cast<geometric::PathGeometric&>(inpath).getStates()
		: projectedWaypoints;
	unsigned int numSegments = numStates - 1;
    OMPL_WARN("Projecting a path with %d waypoints onto constraint manifold", numStates);
	for (unsigned int i = 0; i < numSegments; ++i)
	{
		outpath.append(cloneState(instates[i]));
		if (!constrainedExtend(instates[i], instates[i + 1], outpath.getStates()))
		{
            OMPL_WARN("Projection of segment %d failed!", i);
			result = false;
			break;
		}
	}
	if (result)
    {
		outpath.append(cloneState(instates[numSegments]));
        OMPL_WARN("Path projection was successful!");

        for (unsigned int i = 0; i < outpath.getStateCount(); ++i)
            assert(isValid(outpath.getState(i)));
    }
	else
	{
		//for (unsigned int i = 0; i < projectedWaypoints.size(); ++i)
		//	freeState(projectedWaypoints[i]);
		// for (unsigned int i = 0; i < outpath.getStateCount(); ++i)
		// 	freeState(outpath.getState(i));
	}

	return result;
}

bool ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath(
    const geometric::PathGeometric &inpath, geometric::PathGeometric &outpath, bool shorten) const
{
    // remove any states in the outpath if necessary
    std::vector<State*> &outstates(outpath.getStates());
    for (unsigned int i = 0; i < outstates.size(); ++i)
        freeState(outstates[i]);
    outstates.clear();

    // check if start of path needs to be fixed
    unsigned int validStart = 0;
    if (!isValid(inpath.getState(0)))
    {
        State* scratchState = cloneState(inpath.getState(0));
        if (ci_->project(scratchState) && isValid(scratchState))
        {
            OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath: fixed up start state");
            outpath.append(inpath.getState(0));
            outpath.append(scratchState);
            freeState(scratchState);
            validStart = 1;
        }
        else
        {
            OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath: failed to fix up start state");
            freeState(scratchState);
            outpath = inpath;
            return false;
        }
    }
    else
        // start state was valid, so add it to outpath
        outpath.append(inpath.getState(0));

    // check if end of path needs to be fixed
    unsigned int numStates = inpath.getStateCount();
    State* endState = const_cast<State*>(inpath.getState(numStates - 1));
    if (!isValid(endState))
    {
        endState = cloneState(endState);
        if (ci_->project(endState) && isValid(endState))
            OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath: fixed up end state");
        else
        {
            OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath: failed to fix up end state");
            freeState(endState);
            outpath = inpath;
            return false;
        }
    }

    bool success = true;
    for (unsigned int i = 1; i < numStates - 1; ++i)
        // project path between pairs of valid states; exit loop when projection fails
        if (isValid(inpath.getState(i)) && !subdivideAndProject(outpath, inpath.getState(i)))
        {
            success = false;
            break;
        }
    if (success)
    {
        // Last state can be fixed-up state or valid inpath end state.
        // In both cases we can skip the isValid call.
        if (!subdivideAndProject(outpath, endState))
            success = false;
        else
            outpath.append(endState);
    }

    if (!success)
    {
        OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath: projection failed");
        outpath = inpath;
        if (endState != inpath.getState(numStates - 1))
            freeState(endState);
        return false;
    }

    // check if we can reduce the number of waypoints in the outpath
    if (shorten)
    {
        unsigned int outNumStates = outpath.getStateCount();
        if (outNumStates > 2 + validStart)
        {
            std::vector<State*> &outstates(outpath.getStates()), scratchStates;
            unsigned int from = validStart, to = from + 2;
            const MotionValidatorPtr &mv(getMotionValidator());
            for (unsigned int i = 0; i <= validStart; ++i)
                scratchStates.push_back(cloneState(outstates[i]));
            while (true)
            {
                // greedy approach: find the largest j s.t. the motion i->j is valid
                while (to < outstates.size() && mv->checkMotion(outstates[from], outstates[to])) ++to;
                // this shouldn't happen: the last state in outpath is guaranteed to be valid
                if (to == outstates.size() - 1)
                {
                    OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath: this cannot happen! Let's pretend it didn't.");
                    scratchStates.push_back(cloneState(outstates[to]));
                    break;
                }
                scratchStates.push_back(cloneState(outstates[to - 1]));
                if (to == outstates.size())
                    break;
                from = to - 1;
                to++;
            }
            if (scratchStates.size() < outstates.size())
                // we were able to shorten the path
                std::swap(scratchStates, outstates);
            // free memory
            for (unsigned int i = 0; i < scratchStates.size(); ++i)
                freeState(scratchStates[i]);
        }
    }

    // add last segment (if necessary) between fixed-up end state and inpath endstate.
    if (endState != inpath.getState(numStates - 1))
        outpath.append(inpath.getState(numStates - 1));

    OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProjectPath:");
    OMPL_WARN("\torig. path: %d states, projected path: %d states",
        numStates, outpath.getStateCount());

    return true;
}

bool ompl::base::ConstrainedSpaceInformation::subdivideAndProject(
    geometric::PathGeometric &outpath, const State* waypoint, unsigned int waypointLimit) const
{
    StateSpacePtr ss = getStateSpace();
    if (waypointLimit == 0)
        // put a limit on how much longer the projected path can be
        waypointLimit = outpath.getStateCount() + 10 * ss->validSegmentCount(
            outpath.getStates().back(), waypoint);
    if (outpath.getStateCount() > waypointLimit)
    {
        OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProject: projected path is too long");
        return false;
    }

    static const unsigned int MAX_ATTEMPTS = 20;

    const State *from = outpath.getStates().back(), *to = waypoint;
    if (ss->validSegmentCount(from, to) > 1)
    {
        State* scratchState = allocState();
        bool proj, valid, success = false;
        double dist_midpt = -1, dist_endpt2 = -1;
        for (unsigned int i = 0; i < MAX_ATTEMPTS; ++i)
        {
            ss->interpolate(from, to, 0.5, scratchState);
            if ((proj = ci_->project(scratchState)) &&
                (valid = isValid(scratchState)) &&
                (dist_midpt = distance(from, scratchState)) <= (dist_endpt2 = 2.0*distance(from, to)))
            {
                success = true;
                break;
            }
        }
        if (!success)
        {
            if (!proj)
                OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProject: projection failed");
            else if (!valid)
                OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProject: projected state is not valid");
            else if (dist_midpt > dist_endpt2)
                OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProject: projected state moved too far away: %g > %g",
                    dist_midpt, dist_endpt2);
            else
                OMPL_WARN("ompl::base::ConstrainedSpaceInformation::subdivideAndProject: unknown error");
        }
        else if (subdivideAndProject(outpath, scratchState, waypointLimit - 1))
        {
            outpath.append(scratchState);
            if (subdivideAndProject(outpath, waypoint, waypointLimit))
                outpath.append(waypoint);
            else
                success = false;
        }
        else
            success = false;
         freeState(scratchState);

         return success;
    }
    return true;
}
