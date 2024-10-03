#pragma once

namespace Btwxt {

    // A mutex that just creates a new mutex when you try to copy it. This just allows having a unique mutex class member for each copy of a class without having to declare a copy constructor for that class.
    class copiable_mutex_member : public std::mutex
    {
    public:
	    copiable_mutex_member() = default;
	    copiable_mutex_member([[maybe_unused]] copiable_mutex_member const& other);
	    copiable_mutex_member& operator=([[maybe_unused]] copiable_mutex_member const& other);
    };

}