from tensorGlow.core.index import Index, IndexType

ST = IndexType(name="spacetime", dimension=4, dummy_index=None)
m  = Index("m", ST, up=True)
m_ = Index("m", ST, up=False)  # lowered


class TestIndexType:
    def test_attributes(self):
        assert ST.name == "spacetime"
        assert ST.dimension == 4


class TestIndex:
    def test_default_contravariant(self):
        assert Index("m", ST).up is True

    def test_attributes(self):
        assert m.name == "m"
        assert m.index_type is ST


class TestIndexNegation:
    def test_flips_variance(self):
        assert (-m).up is False
        assert (-m_).up is True

    def test_double_negation_identity(self):
        assert (--m).up is True

    def test_preserves_name_and_type(self):
        assert (-m).name == "m"
        assert (-m).index_type is ST

    def test_returns_new_object(self):
        assert -m is not m
