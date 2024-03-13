## Creating a domain from a generator

Domain(gen::Base.Generator) = generator_domain(gen)

generator_domain(gen::Base.Generator) = generator_domain(gen.f, gen.iter)

# Example: Domain(x for x in 0..1)
generator_domain(f::typeof(identity), iter) = checkdomain(iter)
# Example: Domain(x>1 for x in 0..2)
generator_domain(f, iter) = BoundedIndicatorFunction(f, checkdomain(iter))

# Example: Domain(x for x in 0..2 if x > 1)
generator_domain(f::typeof(identity), iter::Base.Iterators.Filter) =
    generator_domain(iter.flt, iter.itr)
# Example: Domain(x>1 for x in 0..2 if x < 1.5)
generator_domain(f, iter::Base.Iterators.Filter) =
    generator_domain(t -> f(t) && iter.flt(t), iter.itr)
# Example: Domain(x*y>0 for (x,y) in UnitDisk())
generator_domain(f, iter::Base.Iterators.ProductIterator) =
    generator_productdomain(f, iter.iterators)
# avoids ambiguity warning by Aqua, method is never called
generator_domain(f::typeof(identity), iter::Base.Iterators.ProductIterator) =
    generator_productdomain(f, iter.iterators)

function generator_productdomain(f, iterators)
    domain = TupleProductDomain(iterators)
    BoundedIndicatorFunction(f, domain)
end
