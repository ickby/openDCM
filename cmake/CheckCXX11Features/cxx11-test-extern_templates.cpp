
template< typename T >
struct foo {};

extern template class foo<int>;

int main()
{
    return 0;
}

