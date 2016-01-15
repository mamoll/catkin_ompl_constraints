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

	if (!waypointsValid)
	{
		base::State* scratchState = allocState();
		for (unsigned int i = 0; i < numStates; ++i)
		{
			copyState(scratchState, inpath.getState(i));
			if (!ci_->project(scratchState) || !isValid(scratchState))
			{
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
	for (unsigned int i = 0; i < numSegments; ++i)
	{
		outpath.append(cloneState(instates[i]));
		if (!constrainedExtend(instates[i], instates[i + 1], outpath.getStates()))
		{
			result = false;
			break;
		}
	}
	if (result)
		outpath.append(cloneState(instates[numSegments]));
	else
	{
		for (unsigned int i = 0; i < projectedWaypoints.size(); ++i)
			freeState(projectedWaypoints[i]);
		for (unsigned int i = 0; i < outpath.getStateCount(); ++i)
			freeState(outpath.getState(i));
	}

	return result;
}
