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

#ifndef OMPL_BASE_CONSTRAINED_SPACE_INFORMATION_
#define OMPL_BASE_CONSTRAINED_SPACE_INFORMATION_

#include "ompl/base/SpaceInformation.h"
#include "ompl/base/ConstraintInformation.h"
#include "ompl/base/PlannerTerminationCondition.h"
#include <boost/container/slist.hpp>

namespace ompl
{
    namespace geometric
    {
        class PathGeometric;
    }

    namespace base
    {
        /// @cond IGNORE
        /** \brief Forward declaration of ompl::base::ConstrainedSpaceInformation */
        OMPL_CLASS_FORWARD(ConstrainedSpaceInformation);
        /// @endcond

        class ConstrainedSpaceInformation : public SpaceInformation
        {
        public:

            ConstrainedSpaceInformation(const StateSpacePtr &space);

            virtual ~ConstrainedSpaceInformation(void)
            {
            }

            void setConstraintInformation(const ConstraintInformationPtr& ci);

            const ConstraintInformationPtr& getConstraintInformation() const;

            /// \brief Starting at a, interpolate toward b along the constraint manifold.
            /// Store the resulting states in result.
            bool constrainedExtend(const State *a, const State *b,
                                   std::vector<State*>& result) const;

            bool shortcutPath(geometric::PathGeometric *path, const PlannerTerminationCondition &ptc) const;

            bool projectPath(const geometric::PathGeometric &inpath, geometric::PathGeometric &outpath,
                bool waypointsValid = true) const;

            bool subdivideAndProjectPath(const geometric::PathGeometric &inpath, geometric::PathGeometric &outpath) const;

        protected:
            bool subdivideAndProject(boost::container::slist<base::State*> &outpath, const boost::container::slist<base::State*>::iterator &pos) const;

            ConstraintInformationPtr ci_;

            /** \brief The random number generator used by shortcutPath */
            mutable RNG rng_;
        };

    }

}

#endif
